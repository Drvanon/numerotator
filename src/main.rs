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

struct IMGTConservedAminoAcids {
    first_cys: usize,
    conserved_trp: usize,
    hydrophobic_89: usize,
    second_cys: usize,
    j_trp_or_phe: usize,
}

fn count_gaps_in_sequence_before_index(sequence: &[u8], index: usize) -> usize {
    sequence
        .into_iter()
        .take(index)
        .filter(|char| **char == b'-')
        .count()
}

#[derive(Debug)]
enum IMGTError {
    InvalidAlignment,
}

impl IMGTConservedAminoAcids {
    fn is_valid_alignment(alignment: &[u8]) -> bool {
        let (&aa_23, &aa_41, &aa_89, &aa_104, &aa_118) = match alignment
            .into_iter()
            .enumerate()
            .filter_map(|(position, char)| {
                [23, 41, 89, 104, 118]
                    .contains(&(position + 1))
                    .then(|| char)
            })
            .collect_tuple()
        {
            Some(tup) => tup,
            None => return false,
        };

        aa_23 == b'W'
            && aa_41 == b'C'
            && aa_104 == b'C'
            && [b'F', b'W'].contains(&aa_118)
            && [b'A', b'I', b'L', b'M', b'F', b'W', b'Y', b'V'].contains(&aa_89)
    }

    fn from_alignment(alignment: &[u8]) -> Result<Self, IMGTError> {
        if !Self::is_valid_alignment(alignment) {
            return Err(IMGTError::InvalidAlignment);
        }

        Ok(Self {
            first_cys: 23 - count_gaps_in_sequence_before_index(alignment, 23),
            conserved_trp: 41 - count_gaps_in_sequence_before_index(alignment, 41),
            hydrophobic_89: 89 - count_gaps_in_sequence_before_index(alignment, 89),
            second_cys: 104 - count_gaps_in_sequence_before_index(alignment, 104),
            j_trp_or_phe: 118 - count_gaps_in_sequence_before_index(alignment, 118),
        })
    }
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

// TODO: Write a proper stockholm reader.
fn initialize_imgt_alignments() -> HashMap<&'static str, IMGTConservedAminoAcids> {
    let stockholm_data = include_str!("data/reference.stockholm");
    stockholm_data
        .split_ascii_whitespace()
        .tuples()
        .filter_map(|(id, alignment)| {
            IMGTConservedAminoAcids::is_valid_alignment(alignment.as_bytes()).then(|| {
                (
                    id,
                    IMGTConservedAminoAcids::from_alignment(alignment.as_bytes())
                        .expect("Invalid alignment in reference alignments."),
                )
            })
        })
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

    let ref_seqs = initialize_ref_seqs();
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
}

fn initialize_ref_seqs() -> Vec<fasta::Record> {
    fasta::Reader::new(io::Cursor::new(include_bytes!("data/reference.fasta")))
        .records()
        .map(|record_result| record_result.expect("Reference records should be valid."))
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_validity_of_references() {
        let ref_seq = initialize_ref_seqs();
        ref_seq
            .into_iter()
            .for_each(|rec| assert!(IMGTConservedAminoAcids::is_valid_alignment(rec.seq())))
    }

    #[test]
    fn test_imgtconserved_amino_acids_from_str() {
        let test_str = "QVQLVQSGA-EVKKPGASVKVSCKASGYTF----TSYGISWVRQAPGQGLEWMGWISAY--NGNTNYAQKLQ-GRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAR--------MDVWGQGTTVTVSS";
        let conserved_aas = IMGTConservedAminoAcids::from_alignment(test_str.as_bytes()).unwrap();
        assert_eq!(conserved_aas.first_cys, 22);
        assert_eq!(conserved_aas.conserved_trp, 36);
        assert_eq!(conserved_aas.hydrophobic_89, 81);
        assert_eq!(conserved_aas.second_cys, 96);
        assert_eq!(conserved_aas.j_trp_or_phe, 102);
    }
}
