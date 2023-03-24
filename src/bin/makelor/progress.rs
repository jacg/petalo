/// Progress bar and statistics for `makelor` excutable
pub (super) struct Progress(Mutex<Inner>);

struct Inner {
    n_files_given: u16,
    n_files_read: u16,
    n_events_read: u64,
    n_events_processed: u64,
    n_lors_made: u64,
    files_bar: ProgressBar,
    failed_files: Vec<hdf5::Error>,
}

impl Inner {
    fn update(&mut self) {
        let Inner {
            n_files_given,
            n_files_read,
            n_events_read,
            n_events_processed,
            n_lors_made,
            files_bar,
            failed_files,
        } = &self;
        let n_files_failed = failed_files.len();
        files_bar.tick();
        //println!("Files read: {n_files_read}/{n_files_given} (failed: {n_files_failed})");
    }
}

impl Progress {

    pub (super) fn new(infiles: &[PathBuf]) -> Self {
        let bar = ProgressBar::new(infiles.len() as u64).with_message(infiles[0].display().to_string());
        bar.set_style(ProgressStyle::default_bar()
                      .template("[{elapsed_precise}] {wide_bar} {pos}/{len} ({eta_precise})")
                      .unwrap()
        );
        bar.tick();
        Self (
            Mutex::new(
                Inner {
                    n_files_given: infiles.len() as u16,
                    n_files_read: 0,
                    n_events_read: 0,
                    n_events_processed: 0,
                    n_lors_made: 0,
                    files_bar: bar,
                    failed_files: vec![],
                }
            )
        )
    }

    pub (super) fn read_file_start(&self, file: &Path) {
        // Maybe only need read_file_done
    }

    pub (super) fn read_file_done<T>(&self, result: &Result<Vec<T>>) {
        let mut data = self.0.lock().unwrap();
        data.n_files_read += 1;
        if let Err(e) = result {
            data.failed_files.push(e.clone()) // TODO: this should be a 'progress' bar
        }
        data.files_bar.inc(1);
        data.update();
    }

    pub (super) fn grouped<T>(&self, groups: &[T]) {
        let mut data = self.0.lock().unwrap();
        data.n_events_read += groups.len() as u64;
    }

    // pub (super) fn file_succeeded(&self, batch: &LorBatch) {
    //     let mut data = self.0.lock().unwrap();
    //     data.n_lors_made += batch.lors.len() as u64;
    //     todo!()
    // }

    pub (super) fn final_report(&self, lors: &Vec<petalo::io::hdf5::Hdf5Lor>) {
        let data = self.0.lock().unwrap();
        let n_lors = lors.len();
        let n_events = data.n_events_read;
        let percent = (100 * n_lors) as f32 / n_events as f32;
        println!("Found {} LORs in {} events ({percent:.1}%)", group_digits(n_lors), group_digits(n_events));
        if !data.failed_files.is_empty() {
            println!("Warning: the following errors were encountered when reading files input:");
            for err in data.failed_files.iter() {
                println!("  {err}");
            }
            let n = data.failed_files.len();
            let plural = if n == 1 { "" } else { "s" };
            println!("Warning: failed to read {} file{}:", n, plural);
        }
        // self.files_bar.finish_with_message("<finished processing files>");
        // println!("{} / {} ({}%) events produced LORs", group_digits(lors.len()), group_digits(n_events),
        //          100 * lors.len() / n_events);
    }

}

// ----- Imports -----------------------------------------------------------------------------------------
use std::{
    sync::Mutex,
    path::{Path, PathBuf},
};
use hdf5::Result;
use indicatif::{ProgressBar, ProgressStyle};
use petalo::utils::group_digits;
