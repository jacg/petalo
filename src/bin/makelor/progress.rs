/// Progress bar and statistics for `makelor` excutable
pub (super) struct Progress(Mutex<Inner>);

struct Inner {
    n_files_read: u16,
    n_events_read: u64,
    n_lors_made: u64,
    files_bar: ProgressBar,
    failed_files: Vec<hdf5::Error>,
}

impl Inner {
    fn update(&mut self) {
        let Inner { n_events_read, n_lors_made, files_bar, .. } = &self;
        let percent = (100 * n_lors_made) as f32 / *n_events_read as f32;
        files_bar.set_message(format!("Found {:>12} LORs\n in   {:>12} events ({percent:.1}%)",
                                      group_digits(n_lors_made), group_digits(n_events_read)));
        //files_bar.tick();
    }
}

impl Progress {

    pub (super) fn new(infiles: &[PathBuf]) -> Self {
        let bar = ProgressBar::new(infiles.len() as u64).with_message("Making LORs ...");
        bar.set_style(ProgressStyle::default_bar()
                      .template("{msg}\n[{elapsed_precise}] {wide_bar} {pos}/{len} ({eta_precise})")
                      .unwrap()
        );
        bar.tick();
        Self (
            Mutex::new(
                Inner {
                    n_files_read: 0,
                    n_events_read: 0,
                    n_lors_made: 0,
                    files_bar: bar,
                    failed_files: vec![],
                }
            )
        )
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

    pub (super) fn lor(&self, _: &Hdf5Lor) {
        let mut data = self.0.lock().unwrap();
        data.n_lors_made += 1;
    }

    pub (super) fn final_report(&self) {
        let mut data = self.0.lock().unwrap();
        data.update();
        data.files_bar.finish();
        // self.files_bar.finish_with_message("<finished processing files>");
        if !data.failed_files.is_empty() {
            println!("Warning: the following errors were encountered when reading files input:");
            for err in data.failed_files.iter() {
                println!("  {err}");
            }
            let n = data.failed_files.len();
            let plural = if n == 1 { "" } else { "s" };
            println!("Warning: failed to read {} file{}:", n, plural);
        }
    }

}

// ----- Imports -----------------------------------------------------------------------------------------
use std::{
    sync::Mutex,
    path::PathBuf,
};
use hdf5::Result;
use indicatif::{ProgressBar, ProgressStyle};
use petalo::{utils::group_digits, io::hdf5::Hdf5Lor};
