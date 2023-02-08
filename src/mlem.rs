use ndarray::azip;

use rayon::prelude::*;

use units::{
    Length, PerLength, AreaPerMass,
    todo::{Lengthf32, Intensityf32},
    ratio_, mm, kg,
};

use crate::{
    LOR,
    config::mlem::Tof,
    fov::FOV,
    gauss::make_gauss_option,
    image::{Image, ImageData},
    system_matrix::{
        FovHit, FoldState, SystemMatrixRow, Siddon,
        lor_fov_hit, system_matrix_elements,
        project_one_lor, back_project, forward_project,
    },
};

pub static mut N_MLEM_THREADS: usize = 1;

pub fn mlem(fov          : FOV,
            measured_lors: &[LOR],
            tof          : Option<Tof>,
            sensitivity  : Option<Image>,
            n_subsets    : usize,
) -> impl Iterator<Item = (Image, usize, usize)> + '_ {

    // Start off with a uniform image
    let mut image = Image::ones(fov);

    let sensitivity = sensitivity.or_else(|| Some(Image::ones(fov))).unwrap();

    let len = measured_lors.len();
    let set_size = len / n_subsets; // TODO: remainder LORs ignored
    let (mut iteration, mut subset) = (1, 1);

    // Return an iterator which generates an infinite sequence of images,
    // each one made by performing one MLEM iteration on the previous one
    std::iter::from_fn(move || {
        let lo = (subset - 1) * set_size;
        let hi = lo + set_size;
        let (old_iteration, old_subset) = (iteration, subset);
        subset += 1;
        if subset > n_subsets {
            subset = 1;
            iteration += 1;
        }
        one_iteration(&mut image, &measured_lors[lo..hi], &sensitivity.data, tof);
        Some((image.clone(), old_iteration, old_subset)) // TODO see if we can sensibly avoid cloning
    })
}

fn one_iteration(image: &mut Image, measured_lors: &[LOR], sensitivity: &[Intensityf32], tof: Option<Tof>) {

    // -------- Prepare state required by serial/parallel fold --------------

    // TOF adjustment to apply to the weights
    let tof: Option<_> = make_gauss_option(tof);

    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let shared_image = &*image;
    let initial_thread_state = || {
        let (backprojection, system_matrix_row) = projection_buffers(image.fov);
        (backprojection, system_matrix_row, shared_image, &tof)
    };

    // -------- Project all LORs forwards and backwards ---------------------
    let n_mlem_threads = unsafe {
        // SAFETY: modified only once, at the beginning of bin/mlem.rs::main()
        N_MLEM_THREADS
    };
    let job_size = measured_lors.len() / n_mlem_threads;
    let fold_result = measured_lors
        .par_iter()
        // Rayon is too eager in spawning small jobs, each of which requires the
        // construction and subsequent combination of expensive accumulators
        // (whole `Image`s). So here we try to limit it to one job per thread.
        .with_min_len(job_size)
        .with_max_len(job_size)
        .fold(initial_thread_state, project_one_lor);

    // -------- extract relevant information (backprojection) ---------------
    let backprojection = fold_result
        // Keep only the backprojection (ignore weights and indices)
        .map(|tuple| tuple.0)
        // Sum the backprojections calculated on each thread
        .reduce(|| Image::zeros_buffer(image.fov), elementwise_add);

    // -------- Correct for attenuation and detector sensitivity ------------
    apply_sensitivity_image(&mut image.data, &backprojection, sensitivity);
}

pub fn projection_buffers(fov: FOV) -> (ImageData, SystemMatrixRow) {
    // Allocating these anew for each LOR had a noticeable runtime cost, so we
    // create them up-front and reuse them.
    (
        // The ackprojection (or sensitivity image) being constructed in a given
        // MLEM current_iteration (or sensitivity image calculation)
        Image::zeros_buffer(fov),

        // Weights and indices are sparse storage of the slice through the
        // system matrix which corresponds to the current
        Siddon::buffer(fov)
    )
}

fn elementwise_add(a: Vec<f32>, b: Vec<f32>) -> Vec<f32> {
    a.iter().zip(b.iter()).map(|(l,r)| l+r).collect()
}

fn sensitivity_one_lor<'i, 'g, G>(state: FoldState<'i, 'g, G>, lor: LOR) -> FoldState<'i, 'g, G>
where
    G: Fn(Length) -> PerLength
{
    let (mut backprojection, mut system_matrix_row, attenuation, tof) = state;

    // Need to return the state from various match arms
    macro_rules! return_state { () => (return (backprojection, system_matrix_row, attenuation, tof)); }

    // Find active voxels (slice of system matrix) WITHOUT TOF
    // Analyse point where LOR hits FOV
    match lor_fov_hit(&lor, attenuation.fov) {

        // LOR missed FOV: nothing to be done
        None => return_state!(),

        // Data needed by `system_matrix_elements`
        Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

            // Throw away previous LOR's values
            system_matrix_row.clear();

            // Find active voxels and their weights
            system_matrix_elements(
                &mut system_matrix_row,
                next_boundary, voxel_size,
                index, delta_index, remaining,
                tof_peak, tof
            );

            // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
            for (i, _) in &system_matrix_row {
                if i >= backprojection.len() { return_state!(); }
            }

            let integral = forward_project(&system_matrix_row, attenuation);
            let attenuation_factor = (-integral).exp();
            // Backprojection of LOR onto sensitivity image
            back_project(&mut backprojection, &system_matrix_row, attenuation_factor);
            return_state!();
        }
    }
}

// Too much copy-paste code reuse from project_one_lor. This is because the
// latter (and the functions it uses) was heavily optimized, at the cost of ease
// of reuse.

/// Create sensitivity image by backprojecting LORs. In theory this should use
/// *all* possible LORs. In practice use a representative sample.
pub fn sensitivity_image(density: Image, lors: impl ParallelIterator<Item = LOR>, n_lors: usize, rho_to_mu: AreaPerMass) -> Image {
    // Convert from [density in kg/m^3] to [mu in mm^-1]
    let rho_to_mu: f32 = ratio_({
        let kg = kg(1.0);
        let  m = mm(1000.0);
        let rho_unit = kg / (m * m * m);
        let  mu_unit = 1.0 / mm(1.0);
        rho_to_mu / (mu_unit / rho_unit)
    });
    let mut attenuation = density;
    for voxel in &mut attenuation.data {
        *voxel *= rho_to_mu;
    }

    // TOF should not be used as LOR attenuation is independent of decay point
    let notof = make_gauss_option(None);

    // Closure preparing the state needed by `fold`: will be called by
    // `fold` at the start of every thread that is launched.
    let initial_thread_state = || {
        let (backprojection, system_matrix_row) = projection_buffers(attenuation.fov);
        (backprojection, system_matrix_row, &attenuation, &notof)
    };

    // -------- Project all LORs forwards and backwards ---------------------
    let fold_result = lors
        .fold(initial_thread_state, sensitivity_one_lor);

    // -------- extract relevant information (backprojection) ---------------
    let mut backprojection = fold_result
    // Keep only the backprojection (ignore weights and indices)
        .map(|tuple| tuple.0)
    // Sum the backprojections calculated on each thread
        .reduce(|| Image::zeros_buffer(attenuation.fov), elementwise_add);

    // TODO: Just trying an ugly hack for normalizing the image. Do something sensible instead!
    let size = n_lors as f32;
    for e in backprojection.iter_mut() {
        *e /= size
    }
    Image::new(attenuation.fov, backprojection)
}

fn apply_sensitivity_image(image: &mut ImageData, backprojection: &[Lengthf32], sensitivity: &[Intensityf32]) {
    //  TODO express with Option<matrix> and mul reciprocal
    // Apply Sensitivity matrix
    azip!((voxel in image, &b in backprojection, &s in sensitivity) {
        if s > 0.0 { *voxel *= b * s }
        else       { *voxel  = 0.0   }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use units::{mm, mm_, ns, ratio, turn, turn_, Angle, Ratio};
    use rstest::{rstest, fixture};
    use float_eq::assert_float_eq;

    /// Representation of a straght line in 2D
    /// ax + by + c = 0
    #[derive(Debug, Copy, Clone, PartialEq)]
    struct Line { a: Ratio, b: Ratio, c: Length }

    impl Line {
        /// Construct a line passing through `(x,y)` at `angle` to the positive x-axis.
        fn from_point_and_angle((x,y): (Length, Length), angle: Angle) -> Line {
            let (b,a) = inverse_atan2(angle);
            let b = -b;
            Self { a, b, c: -(a*x + b*y) }
        }

        /// The points at which this line crosses a circle of radius `r`,
        /// centred on the origin.
        fn circle_intersection(self, r: Length) -> Points {
            let eps = mm(000.1) * mm(000.1);
            let Self{ a, b, c } = self;
            let a2_b2 = a*a + b*b;
            let x = - (a*c) / a2_b2;
            let y = - (b*c) / a2_b2;
            if  c*c > r*r*a2_b2        + eps { return Points::Zero }
            if (c*c - r*r*a2_b2).abs() < eps { return Points::One { x,  y } }
            let d = r*r - c*c / a2_b2;
            let m = (d / a2_b2).sqrt();
            let (x1, y1) = (x + m*b, y - m*a);
            let (x2, y2) = (x - m*b, y + m*a);
            Points::Two { x1, y1, x2, y2 }
        }
    }

    impl std::fmt::Display for Line {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            let Line { a, b, c } = *self;
            let (a, b, c) = (ratio_(a), ratio_(b), mm_(c));
            write!(f, "{a}x {b:+}y {c:+} = 0")
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq)]
    enum Points { Zero, One{ x: Length, y: Length }, Two{ x1: Length, y1: Length, x2: Length, y2: Length } }

    #[rstest(/**/   x ,   y ,  turns ,      a ,   b ,   c ,
             case( 0.0,  0.0,   0.0  ,     0.0, -1.0,  0.0), // y =  0
             case( 9.0,  0.0,   0.0  ,     0.0, -1.0,  0.0), // y =  0
             case( 0.0,  0.0,   0.25 ,     1.0,  0.0,  0.0), // x =  0
             case( 0.0,  3.0,   0.25 ,     1.0,  0.0,  0.0), // x =  0
             case( 0.0,  2.0,   0.0  ,     0.0, -1.0,  2.0), // y =  2
             case( 5.0,  2.0,   0.0  ,     0.0, -1.0,  2.0), // y =  2
             case( 0.0, -2.0,   0.0  ,     0.0, -1.0, -2.0), // y = -2
             case( 2.0,  0.0,   0.25 ,     1.0,  0.0, -2.0), // x =  2
             case( 2.0,  9.0,   0.25 ,     1.0,  0.0, -2.0), // x =  2
             case(-2.0,  0.0,   0.25 ,     1.0,  0.0,  2.0), // x = -2
             case( 0.0,  0.0,   0.125,     1.0, -1.0,  0.0), // y =  x
             case( 2.0,  0.0,   0.125,     1.0, -1.0, -2.0), // y =  x - 2
             case( 0.0, -2.0,   0.125,     1.0, -1.0, -2.0), // y =  x - 2
             case(-2.0,  0.0,   0.125,     1.0, -1.0,  2.0), // y =  x + 2
             case( 0.0,  2.0,   0.125,     1.0, -1.0,  2.0), // y =  x + 2
             //case( 0.0,  0.0,   0.17618,   1.0, -0.5,  0.0), // y = 2x
    )]
    fn construct_line(x: f32, y: f32, turns: f32, a: f32, b: f32, c: f32) {
        let (x, y, turns) = (mm(x), mm(y), turn(turns));
        let line = Line::from_point_and_angle((x,y), turns);
        let expected = Line { a: ratio(a), b: ratio(b), c: mm(c) };
        println!("constructed: {line}\nexpected   : {expected}");
        assert_float_eq!(ratio_(line.a), a, ulps <= 1);
        assert_float_eq!(ratio_(line.b), b, ulps <= 1);
        assert_float_eq!(   mm_(line.c), c, ulps <= 1);
    }

    #[rstest(/**/   x ,  y , turns,     r   , expected,
             // Vertical lines, close to y = 2
             case( 6.6, 2.1, 0.0  , mm(2.0) , Points::Zero)                                                       ,
             case( 6.6, 2.0, 0.0  , mm(2.0) , Points::One { x : mm( 0.0), y : mm( 2.0)                           }),
             //case( 6.6, 1.9, 0.0  , mm(2.0) , Points::Two { x1: mm(-0.6), y1: mm( 1.9), x2: mm(0.6), y2: mm(1.9) }),
             case( 6.6, 0.0, 0.0  , mm(2.0) , Points::Two { x1: mm(-2.0), y1: mm( 0.0), x2: mm(2.0), y2: mm(0.0) }),
             // Horizontal lines, close to x = 2
             case( 2.1, 9.9, 0.25 , mm(2.0) , Points::Zero)                                                       ,
             case( 2.0, 9.9, 0.25 , mm(2.0) , Points::One { x : mm( 2.0), y : mm( 0.0)                           }),
             //case( 1.9, 9.9, 0.25 , mm(2.0) , Points::Two { x1: mm( 1.9), y1: mm(-0.6), x2: mm(1.9), y2: mm(0.6) }),
             case( 0.0, 9.9, 0.25 , mm(2.0) , Points::Two { x1: mm( 0.0), y1: mm(-2.0), x2: mm(0.0), y2: mm(2.0) }),
    )]
    fn intersect_line_circle(x: f32, y: f32, turns: f32, r: Length, expected: Points) {
        let (x, y, turns) = (mm(x), mm(y), turn(turns));
        let line = Line::from_point_and_angle((x,y), turns);
        let points = line.circle_intersection(r);
        println!("---------> {line} <------------");
        assert_eq!(points, expected);
    }

    /// A version of tan that avoids problems near the poles of tan
    /// angle -> a,b ∈ \[-1, 1\] such that b/a = tan(angle)
    fn inverse_atan2(mut angle: Angle) -> (Ratio, Ratio) {
        let full   : Angle = turn(1.0);
        let half   : Angle = full    / 2.0;
        let quarter: Angle = half    / 2.0;
        let eighth : Angle = quarter / 2.0;
        // Move to preferred range [-1/8, +3/8) turns
        while angle >= 3.0 * eighth {
            angle -= half;
        }
        if angle < eighth { (ratio(1.0)             , angle.tan()) }
        else              { ((quarter - angle).tan(), ratio(1.0) ) }

    }

    /// Tangent of n / 32 full turns
    fn t32(n: i32) -> f32 { ratio_(turn(n as f32 / 32.0).tan()) }

    #[rstest(/**/    turns   ,       x ,       y  ,
             case( 0.0 / 32.0,      1.0,      0.0),
             case( 1.0 / 32.0,      1.0,  t32(1) ),
             case( 2.0 / 32.0,      1.0,  t32(2) ),
             case( 3.0 / 32.0,      1.0,  t32(3) ),
             case( 4.0 / 32.0,      1.0,      1.0),
             case( 5.0 / 32.0,  t32(3) ,      1.0),
             case( 6.0 / 32.0,  t32(2) ,      1.0),
             case( 7.0 / 32.0,  t32(1) ,      1.0),
             case( 8.0 / 32.0,      0.0,      1.0),
             case( 9.0 / 32.0, -t32(1) ,      1.0),
             case(10.0 / 32.0, -t32(2) ,      1.0),
             case(11.0 / 32.0, -t32(3) ,      1.0),
             case(12.0 / 32.0,      1.0,     -1.0),
             case(13.0 / 32.0,      1.0, -t32(3) ),
             case(14.0 / 32.0,      1.0, -t32(2) ),
             case(15.0 / 32.0,      1.0, -t32(1) ),
             case(16.0 / 32.0,      1.0,      0.0),
             case(17.0 / 32.0,      1.0,  t32(1) ),
             case(18.0 / 32.0,      1.0,  t32(2) ),
             case(19.0 / 32.0,      1.0,  t32(3) ),
             case(20.0 / 32.0,      1.0,      1.0),
             case(21.0 / 32.0,  t32(3) ,      1.0),
             case(22.0 / 32.0,  t32(2) ,      1.0),
             case(23.0 / 32.0,  t32(1) ,      1.0),
             case(24.0 / 32.0,      0.0,      1.0),
             case(25.0 / 32.0, -t32(1) ,      1.0),
             case(26.0 / 32.0, -t32(2) ,      1.0),
             case(27.0 / 32.0, -t32(3) ,      1.0),
             case(28.0 / 32.0,      1.0,     -1.0),
             case(29.0 / 32.0,      1.0, -t32(3)),
             case(30.0 / 32.0,      1.0, -t32(2)),
             case(31.0 / 32.0,      1.0, -t32(1)),
    )]
    fn inv_atan2(turns: f32, x: f32, y: f32) {
        let (a,b) = inverse_atan2(turn(turns));
        assert_float_eq!((ratio_(a),ratio_(b)), (x,y), abs <= (3e-7, 4e-7));
    }

    /// `n` uniformly angularly distributed LOR passing through `(x,y)`
    fn n_lors_through(n: usize, (x, y): (Length, Length)) -> Vec<LOR> {
        let mut lors = vec![];
        for angle in n_angles_around_half_circle_starting_from(n, turn(0.01)) {
            let lor = Line::from_point_and_angle((x,y), angle);
            match lor.circle_intersection(DETECTOR_RADIUS) {
                Points::Two { x1, y1, x2, y2 } => {
                    lors.push(LOR::from_components((ns(0.0), ns(0.0)),
                                                   (x1, y1, mm(0.0)),
                                                   (x2, y2, mm(0.0)),
                                                   ratio(1.0)))
                },
                _ => panic!("LOR does not cross detector at two points.")
            }
        }
        lors
    }

    fn n_angles_around_half_circle_starting_from(n: usize, start: Angle) -> impl Iterator<Item = Angle> {
        let mut i = 0;
        let step = turn(0.5) / (n as f32);
        std::iter::from_fn(move || {
            if i < n {
                let angle = start + i as f32 * step;
                i += 1;
                Some(angle)
            } else {
                None
            }
        })
    }

    #[rstest(/**/ n, start    , expected,
             case(0, turn(0.4), vec![]),
             case(1, turn(0.3), vec![turn(0.3)]),
             case(2, turn(0.2), vec![turn(0.2), turn(0.45)]),
             case(2, turn(0.1), vec![turn(0.1), turn(0.35)]),
             case(5, turn(0.2), vec![turn(0.2), turn(0.3), turn(0.4), turn(0.5), turn(0.6)]),
    )]
    fn distributed_angles(n: usize, start: Angle, expected: Vec<Angle>) {
        let angles: Vec<f32> = n_angles_around_half_circle_starting_from(n, start)
            .map(turn_)
            .collect();
        let expected: Vec<f32> = expected
            .into_iter()
            .map(turn_)
            .collect();
        assert_float_eq!(angles, expected, ulps_all <= 1);
    }

    /// Generate points on a square lattice
    fn grid(xs: (i32, i32), ys: (i32, i32)) -> impl Iterator<Item = (i32, i32)> {
        let mut x = xs.0;
        let mut y = ys.0;
        std::iter::from_fn(move || {
            if y > ys.1 { x += 1; y = ys.0 }
            if x > xs.1 { return None }
            y += 1;
            Some((x, (y - 1)))
        })
    }

    // TODO: this should be reused in bin/mlem:main (report_time complicates it)
    /// Return a function which saves images in the given directory
    fn save_each_image_in(directory: String) -> impl FnMut((Image, usize, usize)) {
        use std::path::PathBuf;
        std::fs::create_dir_all(PathBuf::from(&directory)).unwrap();
        move |(image, iteration, subset)| {
            let image_path = PathBuf::from(format!("{directory}/{iteration:02}-{subset:02}.raw"));
            crate::io::raw::Image3D::from(&image).write_to_file(&image_path).unwrap();
        }
    }

    /// The position of a rectangular region of interest in the x-y plane, with
    /// an associated count of decays associated with the region
    struct Roi {
        x: (i32, i32),
        y: (i32, i32),
        activity: usize,
    }

    /// Create a vector of coincidences generated by decays in the foreground
    /// and background ROIs. The background only contributes decays in those
    /// voxels which are not covered by any of the foreground ROIs. To
    /// facilitate generating higher statistics, `scale` is applied to each
    /// activity.
    fn trues_from_rois(foreground_rois: &[&Roi], background_roi: &Roi, scale: usize) -> Vec<LOR> {
        let mut activity = std::collections::HashMap::new();
        // Sum activities in voxels covered by foreground_rois
        for roi in foreground_rois {
            for (x,y) in grid(roi.x, roi.y) {
                *activity.entry((x, y)).or_insert(0) += roi.activity * scale;
            }
        }
        // Set activity in background_roi only in voxels which have not been set by the foregrounds
        for (x,y) in grid(background_roi.x, background_roi.y) {
            activity.entry((x, y)).or_insert(background_roi.activity * scale);
        }
        let mut lors = vec![];
        for ((x,y), a) in activity {
            lors.extend(n_lors_through(a, (mm(x as f32), mm(y as f32))))
        }
        lors
    }

    /// Use a rough model to apply a scattering perturbation to `lors`
    ///
    /// The given LORs are adjusted by randomly picking one end and shifting it
    /// around the detector ring by a random angle within `max_deviation`.
    ///
    /// The number of scatters generated per input LOR is controlled by `scale`.
    fn scatter_lors(lors: impl Iterator<Item = LOR>, max_deviation: Angle, scale: i16) -> Vec<LOR> {
        let max_deviation: f32 = turn_(max_deviation);
        use geometry::Point;
        use rand::prelude::*;

        fn rotate_by(delta: Angle, Point { x, y, z }: Point) -> Point {
            let new_angle = y.atan2(x) + delta;
            let x = DETECTOR_RADIUS * new_angle.cos();
            let y = DETECTOR_RADIUS * new_angle.sin();
            Point { x, y, z }
        }

        // How to scatter one LOR
        let scatter = |mut lor: LOR| {
            let mut rng = thread_rng();
            let delta = turn(rng.gen_range(-max_deviation..=max_deviation));
            if rng.gen() { lor.p1 = rotate_by(delta, lor.p1); }
            else         { lor.p2 = rotate_by(delta, lor.p2); }
            lor
        };

        // How to produce multiple scatters derived single LOR
        let scaled_scatter = |lor: LOR| {
            let mut i = 0;
            std::iter::from_fn(move || {
                if i < scale {
                    i += 1;
                    Some(scatter(lor))
                } else {
                    None
                }
            })
        };

        lors.into_iter()
            .flat_map(scaled_scatter)
            .collect()
    }

    use units::in_base_unit;
    pub const DETECTOR_RADIUS: Length = in_base_unit!(50.0);

    // Regions of interest for re-use in tests
    const ACT_1: usize =   0;
    const ACT_2: usize = 100;
    const ACT_3: usize =  80;
    const BG   : usize =  20;
    #[fixture] fn roi_1() -> Roi { Roi { x: (-15,-10), y: (  9,14), activity: ACT_1 } }
    #[fixture] fn roi_2() -> Roi { Roi { x: ( 10, 15), y: (  6,11), activity: ACT_2 } }
    #[fixture] fn roi_3() -> Roi { Roi { x: (- 7,  7), y: ( -8,-4), activity: ACT_3 } }
    #[fixture] fn roi_b() -> Roi { Roi { x: (-20, 20), y: (-20,20), activity: BG    } }

    #[fixture]
    fn fov() -> FOV {
        let n = 51;
        let l = mm(n as f32);
        FOV::new((l, l, mm(1.0)),
                 (n, n,    1   ))
    }

    #[fixture]
    #[once]
    fn trues_and_scatters(roi_1: Roi, roi_2: Roi, roi_3: Roi, roi_b: Roi) -> (Vec<LOR>, Vec<LOR>) {
        // Generate scatters and trues
        let trues = trues_from_rois(&[&roi_1, &roi_2, &roi_3], &roi_b, 1);
        let noise = scatter_lors(trues.iter().cloned(), turn(0.1), 1);

        // // Write LORs to file, for use in testing of CLI scattergram specification
        // let mut hdf5lors: Vec<Hdf5Lor> = vec![];
        // {
        //     let mklor = | &LOR { p1, p2, .. }, prompt | {
        //         let (E1, E2) = match prompt {
        //             Prompt::True    => (511.0, 511.0),
        //             Prompt::Scatter => (450.0, 450.0),
        //             _ => panic!("Not expecting randoms"),
        //         };
        //         Hdf5Lor {
        //             dt: 0.0,
        //             x1: mm_(p1.x), y1: mm_(p1.y), z1: mm_(p1.z),
        //             x2: mm_(p2.x), y2: mm_(p2.y), z2: mm_(p2.z),
        //             q1: f32::NAN, q2: f32::NAN,
        //             E1, E2
        //         }
        //     };
        //     for lor in &trues { hdf5lors.push(mklor(lor, Prompt::True   )); }
        //     for lor in &noise { hdf5lors.push(mklor(lor, Prompt::Scatter)); }
        // }
        // let filename = std::path::PathBuf::from("test-mlem-images/lors.h5");
        // std::fs::create_dir_all(filename.parent().unwrap()).unwrap();
        // hdf5::File::create(filename).unwrap()
        //     .create_group("reco_info").unwrap()
        //     .new_dataset_builder()
        //     .with_data(&hdf5lors)
        //     .create("lors").unwrap();

        (trues, noise)
    }

    use crate::lorogram::{BuildScattergram as Sc, Prompt};

    #[rstest(/**/ name        , correction,
             case("corr-none" , Sc::new()                                        ),
             case("corr-r"    , Sc::new()             .r_bins(20).r_max(mm(30.0))),
             case("corr-phi"  , Sc::new().phi_bins(20)                           ),
             case("corr-r-phi", Sc::new().phi_bins(20).r_bins(20).r_max(mm(30.0))),
    )]
    fn mlem_reco(correction: Sc, name: &str, fov: FOV,
                 trues_and_scatters: &(Vec<LOR>, Vec<LOR>))
    {
        let (trues, noise) = trues_and_scatters;

        let mut sgram = correction.build();


        // Fill scattergam, if required
        if let Some(sgram) = &mut sgram {
            for lor in trues { sgram.fill(Prompt::True   , lor); }
            for lor in noise { sgram.fill(Prompt::Scatter, lor); }
        }

        // Combines trues and scatters into single collection of prompts
        let mut lors = trues.clone();
        lors.extend(noise);

        // Annotate each LOR with additive correction taken from scattergam
        if let Some(sgram) = sgram {
            for mut lor in &mut lors {
                lor.additive_correction = sgram.value(lor);
            }
        }

        // Perform MLEM reconstruction, saving images to disk
        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        pool.install(|| {
            mlem(fov, &lors, None, None, 1)
                .take(10)
                .for_each(save_each_image_in(format!("test-mlem-images/{name}/")));
        });
        //assert!(false);
    }
}
