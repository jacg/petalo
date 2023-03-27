mod cli;
mod progress;
mod continuous;
mod discrete;

fn main() -> hdf5::Result<()> {
    let args = Cli::parse();

    // Before starting the potentially long computation, make sure that we can
    // write the result to the requested destination. If the directory where
    // results will be written does not exist yet, make it. Panic if that's impossible.
    std::fs::create_dir_all(PathBuf::from(&args.out).parent().unwrap())
        .unwrap_or_else(|e| panic!("\n\nCan't write to \n\n   {}\n\n because \n\n   {e}\n\n", args.out.display()));
    // --- Progress bar --------------------------------------------------------------
    let progress = progress::Progress::new(&args.infiles);
    // --- Process input files -------------------------------------------------------
    let pool = rayon::ThreadPoolBuilder::new().num_threads(args.threads).build().unwrap();
    let xyzs = io::hdf5::sensors::read_sensor_map(&args.infiles[0])?;
    let files = args.infiles.clone();
    macro_rules! go { ($a:expr, $b:expr, $c:expr) => { compose_steps(files, &progress, $a, $b, $c) }; }
    let lors: Vec<_> = pool.install( || match &args.reco {
        d@Reco::Discrete { .. } => go!(read_vertices, group_vertices, lor_from_discretized_vertices(d)),
        Reco::FirstVertex       => go!(read_vertices, group_vertices, lor_from_first_vertices),
        Reco::BaryVertex        => go!(read_vertices, group_vertices, lor_from_barycentre_of_vertices),
        &Reco::Half { q }       => go!(read_qts     , group_qts(q)  , lor_from_hits(&xyzs)),
        &Reco::Dbscan { q, min_count, max_distance } =>
            go!(read_qts, group_qts(q), lor_from_hits_dbscan(&xyzs, min_count, max_distance)),
        Reco::SimpleRec { .. } => todo!(), // Complex process related to obsolete detector design
    });

    // --- Report any files that failed no be read -----------------------------------
    progress.final_report();
    // --- write lors to hdf5 --------------------------------------------------------
    println!("Writing LORs to {}", args.out.display());
    hdf5::File::create(&args.out)?
        .create_group("reco_info")? // TODO rethink all the HDF% group names
        .new_dataset_builder()
        .with_data(&lors)
        .create("lors")?;
    Ok(())
}

fn compose_steps<T, E, G, M>(
    files: Vec<PathBuf>,
    stats: &Progress,
    extract_rows_from_table: E,
    group_by_event         : G,
    make_one_lor           : M,
) -> Vec<Hdf5Lor>
where
    T: Send,
    E: Fn(&Path)  -> hdf5::Result<Vec<T>> + Send + Sync,
    G: Fn(Vec<T>) -> Vec<Vec<T>>          + Send + Sync,
    M: Fn(  &[T]) -> Option<Hdf5Lor>      + Send + Sync,
{
    files
        .into_iter()
        .par_bridge() // With only 2 threads, locking on every LOR is cheap: reconsider if we ever overcome HDF5 global lock
        .map(|file| extract_rows_from_table(&file))      .inspect(|result| stats.read_file_done(result))
        .map(|rows| rows.unwrap_or_else(|_| vec![]))
        .map(group_by_event)                             .inspect(|x| stats.grouped(x))
        .flat_map_iter(|x| x.into_iter()
                       .filter_map(|x| make_one_lor(&x))).inspect(|x| stats.lor(x))
        .collect()
}

fn group_vertices(vertices: Vec<Vertex>) -> Vec<Vec<Vertex>> {
    group_by(|v| v.event_id, vertices)
}

fn group_qts(q: u32) -> impl Fn(Vec<QT>) -> Vec<Vec<QT>> {
    move |qts| {
        qts.into_iter()
           .filter(|h| h.q >= q)
           .group_by(|h| h.event_id)
           .into_iter()
           .map(|(_, group)| group.collect())
           .collect()
    }
}

fn group_by<T>(group_by: impl FnMut(&T) -> u32, qts: impl IntoIterator<Item = T>) -> Vec<Vec<T>> {
    qts.into_iter()
        .group_by(group_by)
        .into_iter()
        .map(|(_, group)| group.collect())
        .collect()
}

fn read_vertices(infile: &Path) -> hdf5::Result<Vec<Vertex>> { io::hdf5::mc::read_vertices(infile, Bounds::none()) }
fn read_qts     (infile: &Path) -> hdf5::Result<Vec<  QT  >> { io::hdf5::sensors::read_qts(infile, Bounds::none()) }

fn vertices_in_scintillator(vertices: &[Vertex]) -> impl Iterator<Item = Vertex> + '_ {
    vertices.iter().cloned().filter(|v| v.volume_id == 0)
}

// ----- Imports -----------------------------------------------------------------------------------------
use std::path::{Path, PathBuf};
use clap::Parser;
use itertools::Itertools;
use rayon::prelude::*;
use petalo::{
    config::mlem::Bounds,
    io::{self,
         hdf5::{
             Hdf5Lor,
             sensors::{QT, SensorMap},
             mc::Vertex
         }
    },
};
use cli::{Cli, Reco};
use progress::Progress;
use continuous::{
    lor_from_first_vertices, lor_from_barycentre_of_vertices, lor_from_hits, lor_from_hits_dbscan,
};
use discrete::lor_from_discretized_vertices;
