use std::error::Error;

use csv;


use serde::Deserialize;
type F = petalo::weights::Length;


use ndarray::prelude::*;

use petalo::weights as pet;

#[derive(Debug, Deserialize)]
struct Event {
    event_id: usize,
    true_r1: F, true_phi1: F, true_z1: F, true_t1: F,
    true_r2: F, true_phi2: F, true_z2: F, true_t2: F,
    reco_r1: F, reco_phi1: F, reco_z1: F,
    reco_r2: F, reco_phi2: F, reco_z2: F,
}

    fn event_to_lor(
        Event{ reco_r1: r1, reco_phi1: phi1, reco_z1: z1, true_t1: t1,
               reco_r2: r2, reco_phi2: phi2, reco_z2: z2, true_t2: t2, .. }: Event)
    -> pet::LOR
{
    let x1 = r1 * phi1.cos();
    let y1 = r1 * phi1.sin();
    let t1 = t1 * 1000.0; // turn into picoseconds

    let x2 = r2 * phi2.cos();
    let y2 = r2 * phi2.sin();
    let t2 = t2 * 1000.0; // turn into picoseconds

    pet::LOR::new(t1, t2,
                  pet::Point::new(x1, y1, z1),
                  pet::Point::new(x2, y2, z2),
    )

}


fn main() -> Result<(), Box<dyn Error>> {

    println!("Float precision: {} bits", petalo::weights::PRECISION);

    use std::time::Instant;

    let mut now = Instant::now();

    let mut report_time = |message| {
        println!("{}: {} ms", message, now.elapsed().as_millis());
        now = Instant::now();
    };

    let mut measured_lors = vec![];

    // Build the CSV reader and iterate over each record.
    let file = std::fs::File::open("run_fastfastmc_1M_events.csv")?;
    let buf_reader = std::io::BufReader::new(file);
    let mut rdr = csv::Reader::from_reader(buf_reader);
    println!("Starting to read events");
    for result in rdr.deserialize() {
        // The iterator yields Result<StringRecord, Error>, so we check the
        // error here..
        let event: Event = result?;
        measured_lors.push(event_to_lor(event));
    }
    report_time("Collected events");

    let vbox = pet::VoxelBox::new((90.0, 90.0, 90.0), (60, 60, 60));

    // TODO: sensitivity matrix, all ones for now
    let sensitivity_matrix = pet::Image::ones(vbox).data;
    // TODO: noise
    let noise = pet::Noise;

    for (n, image) in (pet::Image::mlem(vbox, &measured_lors, &sensitivity_matrix, &noise))
        .take(8)
        .enumerate() {
        report_time("iteration");
        let data: ndarray::Array3<F> = image.data;
        plotit(&data, n).unwrap();
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

    let filename = format!("images/deleteme{:03}.png", n);
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
