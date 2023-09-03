use thiserror::Error;

use std::collections::HashMap;

use bio::{
    alignment::{Alignment, AlignmentOperation},
    io::fasta,
};
use itertools::Itertools;

pub mod annotations;
pub mod regions;
pub mod single_letter_annotations;

/// Container for the positions of a sequence that correspond with IMGT conserved residues in the VREGION.
#[derive(Clone)]
pub struct ConservedResidues {
    first_cys: usize,
    conserved_trp: usize,
    hydrophobic_89: usize,
    second_cys: usize,
    j_trp_or_phe: usize,
}

/// Count the number of gaps in a sequence before a given index.
///
/// Here '-' is the gap character.
pub fn count_gaps_in_sequence_before_index(sequence: &[u8], index: usize) -> usize {
    // TODO: Accept multiple gap characters.
    sequence
        .into_iter()
        .take(index)
        .filter(|char| **char == b'-')
        .count()
}

/// Error for when
#[derive(Debug, Error)]
pub enum IMGTError {
    #[error("Alignment did not have conserved residues in expected places.")]
    InvalidAlignment,

    #[error("Unexpected length ({1}) for region '{0}'.")]
    RegionTooLong(String, usize),

    #[error("CDR3 region too short. Expected at least 5, got {0}")]
    CDR3TooShort(usize),
}

/// Find the position in the alignment sequence that corresponds to a given position.
///
/// For example, say the following alignment exists.
/// ```
///    1 2 3 4 5
/// x: A B C D E
///    | |   | |
/// y: A B - D F
///    1 2   3 4
/// ```
/// Then position 4 in x would correspond to position 3 in y.
fn find_corresponding_position_in_alignment(alignment: &Alignment, pos: usize) -> Option<usize> {
    alignment
        .path()
        .into_iter()
        // TODO: check if there is no off-by-one error here. (Because path starts at one.)
        .find(|(x, _, op)| {
            *x == pos && (*op == AlignmentOperation::Match || *op == AlignmentOperation::Subst)
        })
        .map(|(_, y, _)| y)
}

/// Errors for when transfering conserved residues from one sequence to another.
#[derive(Debug, Error)]
pub enum TransferErr {
    #[error("Conserved residue not in alignment.")]
    ConservedPositionNotInAlignment,
}

impl ConservedResidues {
    /// Confirm that an alignment sequence is a valid IMGT reference sequence.
    ///
    /// The alignment sequence should be of the shap "ABC-DE", similar to what you
    /// would find in a single line of a stockholm file.
    pub fn is_valid_alignment(alignment: &[u8]) -> bool {
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

        aa_23 == b'C'
            && aa_41 == b'W'
            && aa_104 == b'C'
            && [b'F', b'W'].contains(&aa_118)
            && [b'A', b'I', b'L', b'M', b'F', b'W', b'Y', b'V'].contains(&aa_89)
    }

    /// Try to find the conserved residues in a sequence by the alignment.
    ///
    /// Expects that the following relations will be fulfilled:
    ///
    /// - cysteine at position 23
    /// - tryptophan at position 41
    /// - hydrophobic amino acid at position 89
    /// - cysteine at position 104
    /// - phenylalanine or tryptophan at position 118
    pub fn from_alignment(alignment: &[u8]) -> Result<Self, IMGTError> {
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

    /// Identify the conserved residues of a new sequence through the conserved residues of a reference sequence and an alignment between the two.
    pub fn transfer(
        &self,
        alignment: &Alignment,
        _destination: &[u8],
    ) -> Result<Self, TransferErr> {
        // TODO: Ensure that on the destination string, the conserved aas are still there!
        Ok(Self {
            first_cys: find_corresponding_position_in_alignment(alignment, self.first_cys)
                .ok_or(TransferErr::ConservedPositionNotInAlignment)?,
            conserved_trp: find_corresponding_position_in_alignment(alignment, self.conserved_trp)
                .ok_or(TransferErr::ConservedPositionNotInAlignment)?,
            hydrophobic_89: find_corresponding_position_in_alignment(
                alignment,
                self.hydrophobic_89,
            )
            .ok_or(TransferErr::ConservedPositionNotInAlignment)?,
            second_cys: find_corresponding_position_in_alignment(alignment, self.second_cys)
                .ok_or(TransferErr::ConservedPositionNotInAlignment)?,
            j_trp_or_phe: find_corresponding_position_in_alignment(alignment, self.j_trp_or_phe)
                .ok_or(TransferErr::ConservedPositionNotInAlignment)?,
        })
    }
}

/// Load the precomputed and curated vconserved residues associated.
pub fn initialize_conserved_residues() -> HashMap<&'static str, ConservedResidues> {
    // TODO: Write a proper stockholm reader.
    let stockholm_data = include_str!("../data/reference.stockholm");
    stockholm_data
        .split_ascii_whitespace()
        .tuples()
        .filter_map(|(id, alignment)| {
            ConservedResidues::is_valid_alignment(alignment.as_bytes()).then(|| {
                (
                    id,
                    ConservedResidues::from_alignment(alignment.as_bytes())
                        .expect("Invalid alignment in reference alignments."),
                )
            })
        })
        .collect()
}

/// Load the curated reference sequences.
pub fn initialize_ref_seqs() -> Vec<fasta::Record> {
    fasta::Reader::new(std::io::Cursor::new(include_bytes!(
        "../data/reference.fasta"
    )))
    .records()
    .map(|record_result| record_result.expect("Reference records should be valid."))
    .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    use tracing_test::traced_test;
    const TEST_ALIGNMENT_STR: &str = "QVQLVQSGA-EVKKPGASVKVSCKASGYTF----TSYGISWVRQAPGQGLEWMGWISAY--NGNTNYAQKLQ-GRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAR--------MDVWGQGTTVTVSS";

    #[test]
    #[traced_test]
    fn test_valid_conserved_amino_acids() {
        assert!(ConservedResidues::is_valid_alignment(
            TEST_ALIGNMENT_STR.as_bytes()
        ))
    }

    #[test]
    fn test_validity_of_references() {
        let ref_seq = initialize_ref_seqs();
        ref_seq
            .into_iter()
            .for_each(|rec| assert!(ConservedResidues::is_valid_alignment(rec.seq())))
    }

    #[test]
    fn test_imgtconserved_amino_acids_from_str() {
        let conserved_aas =
            ConservedResidues::from_alignment(TEST_ALIGNMENT_STR.as_bytes()).unwrap();
        assert_eq!(conserved_aas.first_cys, 22);
        assert_eq!(conserved_aas.conserved_trp, 36);
        assert_eq!(conserved_aas.hydrophobic_89, 81);
        assert_eq!(conserved_aas.second_cys, 96);
        assert_eq!(conserved_aas.j_trp_or_phe, 102);
    }
}
