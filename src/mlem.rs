use std::path::Path;
use ndarray::azip;

#[cfg(not(feature = "serial"))]
use rayon::prelude::*;

use crate::{io, Lengthf32, Index1_u, Intensityf32};
use crate::{Length, PerLength, Ratio, Time};
use crate::{fov::{lor_fov_hit, FovHit}, system_matrix::{system_matrix_elements, LOR}};
use crate::fov::FOV;
use crate::gauss::make_gauss_option;
use geometry::uom::ratio_;

use crate::image::{Image, ImageData};

impl Image {

    pub fn mlem<'a>(fov: FOV,
                    measured_lors: &'a [LOR],
                    sigma        :     Option<Time>,
                    cutoff       :     Option<Ratio>,
                    sensitivity  :     Option<Self>,
    ) -> impl Iterator<Item = Image> + '_ {

        // Start off with a uniform image
        let mut image = Self::ones(fov);

        let sensitivity = sensitivity.or_else(|| Some(Self::ones(fov))).unwrap();

        // Return an iterator which generates an infinite sequence of images,
        // each one made by performing one MLEM iteration on the previous one
        std::iter::from_fn(move || {
            image.one_iteration(measured_lors, &sensitivity.data, sigma, cutoff);
            Some(image.clone()) // TODO see if we can sensibly avoid cloning
        })
    }

    pub fn from_raw_file(path: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        Ok((&crate::io::raw::Image3D::read_from_file(path)?).into())
    }

    pub fn write_to_raw_file(&self, path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        io::raw::Image3D::from(self).write_to_file(path)?;
        Ok(())
    }

    // Too much copy-paste code reuse from project_one_lor. This is because the
    // latter (and the functions it uses) was heavily optimized, at the cost of
    // ease of reuse.

    // TODO turn this into a method?
    /// Create sensitivity image by backprojecting LORs. In theory this should
    /// use *all* possible LORs. In practice use a representative sample.
    pub fn sensitivity_image(fov: FOV, density: Self, lors: impl Iterator<Item = crate::system_matrix::LOR>, n_lors: usize, stradivarius: f32) -> Self {
        let a = &fov;
        let b = &density.fov;
        if a.n != b.n || a.half_width != b.half_width {
            panic!("For now, attenuation and output image dimensions must match exactly.")
        }
        // TODO convert from density to attenuation coefficient
        let attenuation = density;
        let (mut image, mut weights, mut indices) = projection_buffers(fov);

        // TOF should not be used as LOR attenuation is independent of decay point
        let notof = make_gauss_option(None, None);

        'lor: for lor in lors {
            // Find active voxels (slice of system matrix) WITHOUT TOF
            // Analyse point where LOR hits FOV
            match lor_fov_hit(&lor, fov) {

                // LOR missed FOV: nothing to be done
                None => continue,

                // Data needed by `system_matrix_elements`
                Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

                    // Throw away previous LOR's values
                    weights.clear();
                    indices.clear();

                    // Find active voxels and their weights
                    system_matrix_elements(
                        &mut indices, &mut weights,
                        next_boundary, voxel_size,
                        index, delta_index, remaining,
                        tof_peak, &notof
                    );

                    // Skip problematic LORs TODO: Is the cause more interesting than 'effiing floats'?
                    for i in &indices {
                        if *i >= image.len() { continue 'lor; }
                    }

                    let integral = forward_project(&weights, &indices, &attenuation) / stradivarius;
                    let attenuation_factor = (-integral).exp();
                    //println!("{:<8.2e}  ---  {:<8.2e}", integral, attenuation_factor);
                    // Backprojection of LOR onto image
                    back_project(&mut image, &weights, &indices, attenuation_factor);
                }
            }
        }
        // TODO: Just trying an ugly hack for normalizing the image. Do something sensible instead!
        let size = n_lors as f32;
        for e in image.iter_mut() {
            *e /= size
        }
        Self::new(fov, image)
    }

    fn one_iteration(&mut self, measured_lors: &[LOR], sensitivity: &[Intensityf32], sigma: Option<Time>, cutoff: Option<Ratio>) {

        // -------- Prepare state required by serial/parallel fold --------------

        // TOF adjustment to apply to the weights
        let tof: Option<_> = make_gauss_option(sigma, cutoff);

        // Closure preparing the state needed by `fold`: will be called by
        // `fold` at the start of every thread that is launched.
        let initial_thread_state = || {
            let (backprojection, weights, indices) = projection_buffers(self.fov);
            (backprojection, weights, indices, &self, &tof)
        };

        // Parallel fold takes a function which will return ID value;
        // serial fold takes the ID value itself.
        #[cfg (feature = "serial")]
        // In the serial case, call the function to get one ID value
        let initial_thread_state =  initial_thread_state();

        // Choose between serial parallel iteration
        #[cfg    (feature = "serial") ] let iter = measured_lors.    iter();
        #[cfg(not(feature = "serial"))] let iter = measured_lors.par_iter();

        // -------- Project all LORs forwards and backwards ---------------------

        let fold_result = iter.fold(initial_thread_state, project_one_lor);

        // -------- extract relevant information (backprojection) ---------------

        // In the serial case, there is a single result to unwrap ...
        #[cfg (feature = "serial")]
        let backprojection = fold_result.0; // Keep only backprojection

        // ... in the parallel case, the results from each thread must be
        // combined
        #[cfg(not(feature = "serial"))]
        let backprojection = {
            fold_result
            // Keep only the backprojection (ignore weights and indices)
            .map(|tuple| tuple.0)
            // Sum the backprojections calculated on each thread
            .reduce(|   | zeros_buffer(self.fov),
                    |l,r| l.iter().zip(r.iter()).map(|(l,r)| l+r).collect())
        };

        // -------- Correct for attenuation and detector sensitivity ------------

        apply_sensitivity_image(&mut self.data, &backprojection, sensitivity);

    }

    pub fn ones(fov: FOV) -> Self {
        let [x,y,z] = fov.n;
        let size = x * y * z;
        Self { data: vec![1.0; size], fov}
    }

    pub fn new(fov: FOV, data: ImageData) -> Self {
        let [x, y, z] = fov.n;
        if data.len() != x * y * z {
            // TODO change panic to Option or Result
            panic!("Image data does not match dimensions {:?}", fov.n);
        };
        Image { fov, data }
    }

    pub fn empty(fov: FOV) -> Self {
        let [x,y,z] = fov.n;
        Self::new(fov, vec![0.0; x*y*z])
    }

    pub fn inverted(&self) -> Self {
        let mut inverted = self.clone();
        for e in inverted.data.iter_mut() { *e = 1.0 / *e }
        inverted
    }
}

fn projection_buffers(fov: FOV) -> (ImageData, Vec<Lengthf32>, Vec<usize>) {
    // The backprojection (or sensitivity image) being constructed in a
    // given MLEM iteration (or sensitivity image calculation).
    let image = zeros_buffer(fov);
    // Weights and indices are sparse storage of the slice through the
    // system matrix which corresponds to the current LOR. (Allocating these
    // anew for each LOR had a noticeable runtime cost.)
    let [nx, ny, nz] = fov.n;
    let max_number_of_active_voxels_possible = nx + ny + nz - 2;
    let weights = Vec::with_capacity(max_number_of_active_voxels_possible);
    let indices = Vec::with_capacity(max_number_of_active_voxels_possible);
    (image, weights, indices)
}

// A new empty data store with matching size
fn zeros_buffer(fov: FOV) -> ImageData { let [x,y,z] = fov.n; vec![0.0; x*y*z] }


type FoldState<'r, 'i, 'g, G> = (ImageData , Vec<Lengthf32>, Vec<Index1_u> , &'r &'i mut Image, &'g Option<G>);

fn project_one_lor<'r, 'i, 'g, G>(state: FoldState<'r, 'i, 'g, G>, lor: &LOR) -> FoldState<'r, 'i, 'g, G>
where
    G: Fn(Length) -> PerLength
{
    let (mut backprojection, mut weights, mut indices, image, tof) = state;

    // Analyse point where LOR hits FOV
    match lor_fov_hit(lor, image.fov) {

        // LOR missed FOV: nothing to be done
        None => return (backprojection, weights, indices, image, tof),

        // Data needed by `system_matrix_elements`
        Some(FovHit {next_boundary, voxel_size, index, delta_index, remaining, tof_peak}) => {

            // Throw away previous LOR's values
            weights.clear();
            indices.clear();

            // Find active voxels and their weights
            system_matrix_elements(
                &mut indices, &mut weights,
                next_boundary, voxel_size,
                index, delta_index, remaining,
                tof_peak, tof
            );

            // Forward projection of current image into this LOR
            let projection = forward_project(&weights, &indices, image) * lor.additive_correction;

            // Backprojection of LOR onto image
            back_project(&mut backprojection, &weights, &indices, ratio_(projection));
        }
    }
    // Return updated FoldState
    (backprojection, weights, indices, image, tof)
}

#[inline]
fn forward_project(weights: &[Lengthf32], indices: &[usize], image: &Image) -> Lengthf32 {
    let mut projection = 0.0;
    for (w, &j) in weights.iter().zip(indices.iter()) {
        projection += w * image[j]
    }
    projection
}

#[inline]
fn back_project(backprojection: &mut Vec<Lengthf32>, weights: &[Lengthf32], indices: &[usize], projection: Lengthf32) {
    let projection_reciprocal = 1.0 / projection;
    for (w, &j) in weights.iter().zip(indices.iter()) {
        backprojection[j] += w * projection_reciprocal;
    }
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
    use geometry::{uom::{mm, mm_, ns, ratio, turn, turn_}, Angle};
    use rstest::{rstest, fixture};
    use float_eq::assert_float_eq;

    /// ax + by + c = 0
    #[derive(Debug, Copy, Clone, PartialEq)]
    struct Line { a: Ratio, b: Ratio, c: Length }

    impl Line {
        fn from_point_and_angle((x,y): (Length, Length), angle: Angle) -> Line {
            let (b,a) = inverse_atan2(angle);
            let b = -b;
            let line = Self { a, b, c: -(a*x + b*y) };
            line
        }

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

    /// Add `n` uniformly angularly distributed LOR passing through `(x,y)`, to `lors`
    fn n_decays_at(n: usize, (x, y): (Length, Length), lors: &mut Vec<LOR>) {
        for angle in n_angles_around_half_circle_starting_from(n, turn(0.01)) {
            let lor = Line::from_point_and_angle((x,y), angle);
            match lor.circle_intersection(mm(100.0)) {
                Points::Two { x1, y1, x2, y2 } => {
                    lors.push(LOR::from_components((ns(0.0), ns(0.0)),
                                                   (x1, y1, mm(0.0)),
                                                   (x2, y2, mm(0.0)),
                                                   ratio(1.0)))
                },
                _ => panic!("LOR does not cross detector at two points.")
            }
        }
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

    fn grid(xs: (i32, i32), ys: (i32, i32)) -> impl Iterator<Item = (Length, Length)> {
        let mut x = xs.0;
        let mut y = ys.0;
        return std::iter::from_fn(move || {
            if y > ys.1 { x += 1; y = ys.0 }
            if x > xs.1 { return None }
            y += 1;
            Some((mm(x as f32), mm((y - 1) as f32)))
        })
    }

    // TODO: this should be reused in bin/mlem:main (report_time complicates it)
    fn save_each_image_in(directory: String) -> impl FnMut(&Image) {
        use std::path::PathBuf;
        std::fs::create_dir_all(PathBuf::from(&directory)).unwrap();
        let mut count = 0;
        move |image| {
            count += 1;
            let image_path = PathBuf::from(format!("{directory}/{:02}.raw", count));
            crate::io::raw::Image3D::from(image).write_to_file(&image_path).unwrap();
        }
    }

    struct ROI {
        x: (i32, i32),
        y: (i32, i32),
        activity: usize,
    }

    fn detect_lors(mut lors: &mut Vec<LOR>, rois: &[ROI]) {
        for roi in rois {
            for (x,y) in grid(roi.x, roi.y) { n_decays_at(roi.activity, (x, y), &mut lors) }
        }
    }

    const ACT_1: usize =  60;
    const ACT_2: usize = 100;
    const ACT_3: usize =  80;
    const BG   : usize =  20;
    const NOISE: usize =  30;
    // Regions of interest for re-use in tests
    #[fixture] fn roi_1() -> ROI { ROI { x: (-15,-10), y: (  8,13), activity: ACT_1 - BG } }
    #[fixture] fn roi_2() -> ROI { ROI { x: ( 10, 15), y: (  8,13), activity: ACT_2 - BG } }
    #[fixture] fn roi_3() -> ROI { ROI { x: (- 7,  7), y: ( -8,-4), activity: ACT_3 - BG } }
    #[fixture] fn roi_b() -> ROI { ROI { x: (-20, 20), y: (-20,20), activity: BG         } }
    #[fixture] fn roi_n() -> ROI { ROI { x: (-25, 25), y: (-25,25), activity: NOISE      } }

    #[fixture]
    fn fov() -> FOV {
        let n = 51;
        let l = mm(n as f32);
        FOV::new((l, l, mm(1.0)),
                 (n, n,    1   ))
    }

    enum Bins {
        None,
        R    { nbins    : usize, maxr: Length },
        Phi  { nbins    : usize },
        RPhi { nbins_phi: usize, nbins_r: usize, maxr: Length },
    }

    #[rstest(/**/ name        , bins,
             case("corr-none" , Bins::None),
             case("corr-r"    , Bins::R    {                nbins  : 20, maxr: mm(30.0) }),
             case("corr-phi"  , Bins::Phi  { nbins:     20                              }),
             case("corr-r-phi", Bins::RPhi { nbins_phi: 20, nbins_r: 20, maxr: mm(30.0) }),
    )]
    fn mlem_reco(fov: FOV,
                 roi_1: ROI, roi_2: ROI, roi_3: ROI,
                 roi_b: ROI, roi_n: ROI,
                 bins: Bins, name: &str)
    {

        use crate::lorogram::{Scattergram, Prompt, axis_phi, axis_r};
        use ndhistogram::ndhistogram;

        let mut sgram = match bins {
            Bins::None => None,
            Bins::R { nbins, maxr } =>
                Some(Scattergram::new(&|| Box::new(ndhistogram!(axis_r  (nbins, maxr); usize)))),
            Bins::Phi { nbins } =>
                Some(Scattergram::new(&|| Box::new(ndhistogram!(axis_phi(nbins      ); usize)))),
            Bins::RPhi { nbins_phi, nbins_r, maxr } =>
                Some(Scattergram::new(&|| Box::new(ndhistogram!(
                    axis_phi(nbins_phi      ),
                    axis_r  (nbins_r  , maxr);
                    usize)))),

        };


        let mut trues = vec![]; detect_lors(&mut trues, &[roi_1, roi_2, roi_3, roi_b]);
        let mut noise = vec![]; detect_lors(&mut noise, &[roi_n]);

        if let Some(sgram) = &mut sgram {
            for lor in &trues { sgram.fill(Prompt::True   , lor); }
            for lor in &noise { sgram.fill(Prompt::Scatter, lor); }
        }

        let mut lors = trues;
        lors.extend(noise);

        if let Some(sgram) = sgram {
            for mut lor in &mut lors {
                lor.additive_correction = sgram.value(lor);
            }
        }

        let pool = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();
        let _ = pool.install(|| {
            Image::mlem(fov, &lors, None, None, None)
                .take(10)
                .inspect(save_each_image_in(format!("/tmp/test-mlem-{name}")))
                .for_each(|_| {
                });
        });
        //assert!(false);
    }
}
