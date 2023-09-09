use clap::{value_parser, Parser};
use itertools::Itertools;
use numerotator::imgt::reference::is_valid_alignment;
use tracing::{debug, info, Level};
use tracing_subscriber::FmtSubscriber;

#[derive(Debug, Parser)]
#[command()]
struct Args {
    #[arg(value_parser=value_parser!(std::path::PathBuf))]
    stockholm_file: std::path::PathBuf,

    #[arg(value_parser=value_parser!(std::path::PathBuf))]
    output_fasta_file: std::path::PathBuf,

    #[arg()]
    output_alignments_file: std::path::PathBuf,
}

fn is_alignment_line(line: &str) -> bool {
    (!line.starts_with("//")) && (!line.starts_with("#")) && (!line.is_empty())
}

fn alignment_line_to_fasta(alignment_str: &str) -> String {
    let (id, alignment) = alignment_str
        .split_ascii_whitespace()
        .next_tuple()
        .expect("Alignment string should have id and alignment.");
    let sequence: String = alignment.chars().filter(|char| *char != '-').collect();
    format!(">{}\n{}\n", id, sequence)
}

// TODO: Right now it uses the Anarci
// build_pipeline/curated_alignments/ALL.stockholm file.
// naturaly it should download this itself.
fn main() {
    let args = Args::parse();
    // a builder for `FmtSubscriber`.
    let subscriber = FmtSubscriber::builder()
        // all spans/events with a level higher than TRACE (e.g, debug, info, warn, etc.)
        // will be written to stdout.
        .with_max_level(Level::TRACE)
        // completes the builder.
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");

    info!(
        input_file = args.stockholm_file.as_os_str().to_str(),
        "Reading input file"
    );
    let alignment_data =
        std::fs::read_to_string(args.stockholm_file).expect("Could not open alignments file.");

    debug!(data_size = alignment_data.len(), "Read input file.");

    let mut n_valid_lines = 0;
    // Identify lines with valid sequences.
    let valid_lines: Vec<_> = alignment_data
        .lines()
        .filter(|line| is_alignment_line(line))
        .inspect(|_| n_valid_lines += 1)
        .filter(|line| {
            is_valid_alignment(
                line.split_ascii_whitespace()
                    .skip(1)
                    .take(1)
                    .next()
                    .expect("Valid lines should have at least two parts.")
                    .as_bytes(),
            )
            .is_some()
        })
        .collect();

    debug!(
        n_valid_alignments = valid_lines.len(),
        n_valid_lines, "Validated lines."
    );

    let reference_sequences: String = valid_lines
        .iter()
        .map(|&line| alignment_line_to_fasta(line))
        .collect();

    // Qualifying due to name conflict.
    let reference_alignments: String =
        Itertools::intersperse(valid_lines.into_iter(), "\n").collect();

    std::fs::write(args.output_alignments_file, reference_alignments)
        .expect("Could not write alignments file.");
    std::fs::write(args.output_fasta_file, reference_sequences)
        .expect("Could not write fasta file.");
}
