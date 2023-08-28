use bio::{
    alignment::{self, Alignment},
    io::fasta,
};
use clap::{value_parser, Parser};
use itertools::Itertools;
use numerator::imgt;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command()]
struct Args {
    #[arg(index = 1, num_args=..)]
    sequences: Vec<String>,
    #[arg(short, long, value_parser=value_parser!(PathBuf))]
    sequences_file: Option<PathBuf>,
}

#[derive(Debug)]
pub enum RefSeqErr {
    NoReferenceSequenceFound,
}

struct ReferenceAlignment {
    reference_record: fasta::Record,
    query_record: fasta::Record,
    alignment: Alignment,
}

fn find_best_reference_sequence(
    record: fasta::Record,
    ref_seqs: &Vec<fasta::Record>,
) -> Result<ReferenceAlignment, RefSeqErr> {
    let mut aligner =
        alignment::pairwise::Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });
    ref_seqs
        .into_iter()
        .map(|reference_record| {
            (
                reference_record,
                aligner.local(record.seq(), reference_record.seq()),
            )
        })
        .max_by_key(|(_, alignment)| alignment.score)
        .map(|(reference_record, alignment)| ReferenceAlignment {
            // Cloning here should not be a huge problem, since we only clone once per query sequence.
            reference_record: reference_record.clone(),
            alignment,
            query_record: record,
        })
        .ok_or(RefSeqErr::NoReferenceSequenceFound)
}

fn main() {
    // TODO: make everything u8 based.
    let args = Args::parse();

    let ref_seqs = imgt::initialize_ref_seqs();
    let all_conserved_residues = imgt::initialize_conserved_residues();

    let sequences_from_command_line = args.sequences.into_iter().enumerate().map(|(i, seq)| {
        fasta::Record::with_attrs(
            i.to_string().as_str(),
            Some(format!("sequence {} from the command line", i).as_str()),
            seq[..].as_bytes(),
        )
    });

    let sequences_from_sequence_file = args.sequences_file.and_then(|path| {
        Some(
            fasta::Reader::new(std::fs::File::open(path).expect("Could not open sequences file."))
                .records()
                .map(|record_result| {
                    record_result.expect("Could not parse record in sequences file.")
                }),
        )
    });

    let (reference_alignments, failed_sequences): (Vec<_>, Vec<_>) = sequences_from_command_line
        .chain(sequences_from_sequence_file.into_iter().flatten())
        .map(|query_seq| find_best_reference_sequence(query_seq, &ref_seqs))
        .partition_result();

    let (conserved_sequences, failed_conserved): (Vec<_>, Vec<_>) = reference_alignments
        .into_iter()
        .map(|reference_alignment| {
            all_conserved_residues
                .get(reference_alignment.reference_record.id())
                .expect("Reference sequence id should be in reference alignments aswell.")
                .transfer(
                    &reference_alignment.alignment,
                    reference_alignment.query_record.seq(),
                )
                .and_then(|conserved_sequences| Ok((conserved_sequences, reference_alignment)))
        })
        .partition_result();
}
