use clap::{value_parser, Parser};
use itertools::Itertools;

mod imgt;

#[derive(Parser, Debug)]
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
    !line.starts_with("//") && !line.starts_with("#") && !line.is_empty()
}

fn alignment_line_to_fasta(alignment_str: &str) -> String {
    let (id, alignment) = alignment_str
        .split_ascii_whitespace()
        .next_tuple()
        .expect("Alignment string should have id and alignment.");
    let sequence: String = alignment.chars().filter(|char| *char == '-').collect();
    format!(">{}\n{}\n", id, sequence)
}

// TODO: Right now it uses the Anarci
// build_pipeline/curated_alignments/ALL.stockholm file.
// naturaly it should download this itself.
fn main() {
    let args = Args::parse();
    let alignment_data =
        std::fs::read_to_string(args.stockholm_file).expect("Could not open alignments file.");

    // Identify lines with valid sequences.
    let valid_lines: Vec<_> = alignment_data
        .lines()
        .filter(|line| is_alignment_line(line))
        .filter(|line| imgt::ConservedAminoAcids::is_valid_alignment(line.as_bytes()))
        .collect();

    let reference_sequences: String = valid_lines
        .iter()
        .map(|&line| alignment_line_to_fasta(line))
        .collect();

    let reference_alignments: String = valid_lines.into_iter().collect();

    std::fs::write(args.output_alignments_file, reference_alignments)
        .expect("Could not write alignments file.");
    std::fs::write(args.output_fasta_file, reference_sequences)
        .expect("Could not write fasta file.");
}
