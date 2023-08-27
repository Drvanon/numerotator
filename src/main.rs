use bio::{
    alignment::{self, Alignment},
    io::fasta,
};
use clap::{value_parser, Parser};
use itertools::Itertools;
use std::{collections::HashMap, fs, io, path::PathBuf};

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

fn stockholm_str_to_alignment_ops(stockholm_str: &str) -> Vec<alignment::AlignmentOperation> {
    // x: Current sequence under inspection
    // y: Hypothetical IMGT sequence

    use alignment::AlignmentOperation::*;

    stockholm_str
        .chars()
        .into_iter()
        .map(|char| match char {
            '-' => Del,
            '.' => Del,
            _ => Subst,
        })
        .collect()
}

// TODO: Write a proper stockholm reader.
fn initialize_imgt_alignments() -> HashMap<&'static str, Vec<alignment::AlignmentOperation>> {
    let stockholm_data = include_str!("data/reference.stockholm");
    stockholm_data
        .split_ascii_whitespace()
        .tuples()
        .map(|(id, alignment)| (id, stockholm_str_to_alignment_ops(alignment)))
        .collect()
}

fn alignment_ops_on_x_as_str(operations: &Vec<alignment::AlignmentOperation>, x: &str) -> String {
    let mut cur_index_on_x: usize = 0;
    operations
        .into_iter()
        .map(|operation| match operation {
            alignment::AlignmentOperation::Match => {
                cur_index_on_x += 1;
                x.chars().nth(cur_index_on_x - 1)
            }
            alignment::AlignmentOperation::Subst => {
                cur_index_on_x += 1;
                x.chars().nth(cur_index_on_x - 1)
            }
            alignment::AlignmentOperation::Del => None,
            alignment::AlignmentOperation::Ins => {
                cur_index_on_x += 1;
                x.chars().nth(cur_index_on_x - 1)
            }
            alignment::AlignmentOperation::Xclip(_) => todo!(),
            alignment::AlignmentOperation::Yclip(_) => todo!(),
        })
        .flat_map(std::convert::identity)
        .collect()
}

fn main() {
    // TODO: make everything u8 based.
    let args = Args::parse();

    let ref_seqs: Vec<fasta::Record> =
        fasta::Reader::new(io::Cursor::new(include_bytes!("data/reference.fasta")))
            .records()
            .map(|record_result| record_result.expect("Reference records should be valid."))
            .collect();
    let imgt_alignments = initialize_imgt_alignments();

    let sequences_from_command_line = args.sequences.into_iter().enumerate().map(|(i, seq)| {
        fasta::Record::with_attrs(
            i.to_string().as_str(),
            Some(format!("sequence {} from the command line", i).as_str()),
            seq[..].as_bytes(),
        )
    });

    // TODO: Make this use iterators instead of strings.
    let sequences_from_sequence_file = args.sequences_file.and_then(|path| {
        Some(
            fasta::Reader::new(fs::File::open(path).expect("Could not open sequences file."))
                .records()
                .map(|record_result| {
                    record_result.expect("Could not parse record in sequences file.")
                }),
        )
    });

    let (reference_alignments, failed_sequences): (Vec<_>, Vec<_>) = sequences_from_command_line
        .chain(sequences_from_sequence_file.into_iter().flatten())
        .map(|query_seq| find_best_reference_sequence(query_seq, &ref_seqs))
        .partition(Result::is_ok);

    reference_alignments
        .into_iter()
        .map(Result::unwrap)
        .map(|reference_alignment| {
            let id = reference_alignment.reference_record.id().clone();
            let imgt_operations = imgt_alignments.get(id).expect(
                "Reference sequence in curated sequence should also be in curated alignments.",
            );
            join_alignment_ops(&reference_alignment.alignment.operations, imgt_operations)
        });
}

fn join_alignment_ops(
    from_operations: &Vec<alignment::AlignmentOperation>,
    into_operations: &Vec<alignment::AlignmentOperation>,
) -> Vec<alignment::AlignmentOperation> {
    todo!();
}
