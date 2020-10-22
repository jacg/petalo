// ----------------------------------- CLI -----------------------------------
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "mlem", about = "TODO: describe what this does")]
pub struct Cli {

    /// Number of MLEM iterations to perform
    #[structopt(default_value = "5")]
    n: usize,

}
// --------------------------------------------------------------------------------

use std::error::Error;

use ndarray::prelude::*;

use petalo::weights as pet;

type F = pet::Length;


fn main() -> Result<(), Box<dyn Error>> {

    println!("Float precision: {} bits", petalo::weights::PRECISION);

    use std::time::Instant;

    let mut now = Instant::now();

    let mut report_time = |message| {
        println!("{}: {} ms", message, now.elapsed().as_millis());
        now = Instant::now();
    };

    let filename = "run_fastfastmc_1M_events.bin32";

    let measured_lors = petalo::io::read_lors(filename)?;
    report_time("read bin");

    let vbox = pet::VoxelBox::new((90.0, 90.0, 90.0), (60, 60, 60));

    // TODO: sensitivity matrix, all ones for now
    let sensitivity_matrix = pet::Image::ones(vbox).data;
    // TODO: noise
    let noise = pet::Noise;

    for (n, image) in (pet::Image::mlem(vbox, &measured_lors, &sensitivity_matrix, &noise))
        .take(Cli::from_args().n)
        .enumerate() {
        report_time("iteration");
        let data: ndarray::Array3<F> = image.data;
        plotit(&data, n)?;
        report_time("Plotted");
        // TODO: step_by for print every
    }


    Ok(())
}

use plotters::prelude::*;


fn plotit(image: &ndarray::Array3<F>, n: usize) -> Result<(), Box<dyn std::error::Error>> {

    let slice = image.slice(s![.., .., 29]);
    let size = slice.shape();
    let (nx, ny) = (size[0], size[0]);
    let max = image.fold(0.0, |a:F,b| a.max(*b));

    let filename = format!("images/deleteme{:03}.{}.png", n, pet::PRECISION);
    let root = BitMapBackend::new(&filename, (nx as u32, ny as u32)).into_drawing_area();

    root.fill(&WHITE)?;


    let mut chart = ChartBuilder::on(&root)
        .build_cartesian_2d(0..nx, 0..ny)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    chart.draw_series(
        slice
            .indexed_iter()
            .map(|((x,y), v)| {
                let f = *v / max;
                Rectangle::new(
                    [(x, y), (x + 1, y + 1)],
                    HSLColor(
                        (0.8 - 0.6 * f).into(),
                        0.7,
                        (0.1 + 0.9 * f).into(),
                    )
                        .filled(),
                )
            }),
    )?;

    Ok(())
}
