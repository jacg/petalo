use std::error::Error;
use std::path::PathBuf;
use clap::Parser;
use simset::{compare, custom, standard};

/// Try to parse SimSET (custom) history files
#[derive(Parser, Debug)]
struct Args {
    /// History file to parse
    file: PathBuf,

    /// Second history file for record-by-record comparison
    #[arg(short = 'c', long)]
    compare: Option<PathBuf>,

    /// Maximum number of records to parse
    #[arg(short = 'n', long)]
    stop_after: Option<usize>,

    #[arg(value_enum, short = 't', long = "file-type", default_value_t = HistoryFileType::Custom)]
    history_file_type: HistoryFileType,
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();
    if let Some(ref file2) = args.compare { compare(&args.file, &file2, args.stop_after) }
    else {
        match args.history_file_type {
            HistoryFileType::Standard => standard::show_file(&args.file, args.stop_after),
            HistoryFileType::Custom   =>   custom::show_file(&args.file, args.stop_after),
        }
    }
}

#[derive(clap::ValueEnum, Clone, Debug)]
enum HistoryFileType {
    Standard,
    Custom,
}

impl std::fmt::Display for HistoryFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HistoryFileType::Standard => f.write_str("Standard"),
            HistoryFileType::Custom   => f.write_str("Custom"),
        }
    }
}
