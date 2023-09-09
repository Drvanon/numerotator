use bio::alignment::{Alignment, AlignmentOperation};
use thiserror::Error;

/// Container for the positions of a sequence that correspond with IMGT conserved residues in the VREGION.
#[derive(Clone, Debug)]
pub struct ConservedResidues {
    pub first_cys: usize,
    pub conserved_trp: usize,
    pub hydrophobic_89: usize,
    pub second_cys: usize,
    pub j_trp_or_phe: usize,
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

/// Find the position in the alignment sequence that corresponds to a given position.
///
/// For example, say the following alignment exists.
/// ```txt
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

impl ConservedResidues {
    /// Confirm that an alignment sequence is a valid IMGT reference sequence.
    ///
    /// The alignment sequence should be of the shap "ABC-DE", similar to what you
    /// would find in a single line of a stockholm file.
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

/// Errors for when transfering conserved residues from one sequence to another.
#[derive(Debug, Error)]
pub enum TransferErr {
    #[error("Conserved residue not in alignment.")]
    ConservedPositionNotInAlignment,
}

impl From<&[u8]> for ConservedResidues {
    /// Find the conserved residues in a sequence by the alignment.
    ///
    /// Expects that the following relations will be fulfilled:
    ///
    /// - cysteine at position 23
    /// - tryptophan at position 41
    /// - hydrophobic amino acid at position 89
    /// - cysteine at position 104
    /// - phenylalanine or tryptophan at position 118
    fn from(alignment: &[u8]) -> Self {
        Self {
            first_cys: 23 - count_gaps_in_sequence_before_index(alignment, 23),
            conserved_trp: 41 - count_gaps_in_sequence_before_index(alignment, 41),
            hydrophobic_89: 89 - count_gaps_in_sequence_before_index(alignment, 89),
            second_cys: 104 - count_gaps_in_sequence_before_index(alignment, 104),
            j_trp_or_phe: 118 - count_gaps_in_sequence_before_index(alignment, 118),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    const TEST_ALIGNMENT_STR: &str = "QVQLVQSGA-EVKKPGASVKVSCKASGYTF----TSYGISWVRQAPGQGLEWMGWISAY--NGNTNYAQKLQ-GRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAR--------MDVWGQGTTVTVSS";

    #[test]
    fn test_conserved_amino_acids_from_str() {
        let conserved_aas = ConservedResidues::from(TEST_ALIGNMENT_STR.as_bytes());
        assert_eq!(conserved_aas.first_cys, 22);
        assert_eq!(conserved_aas.conserved_trp, 36);
        assert_eq!(conserved_aas.hydrophobic_89, 81);
        assert_eq!(conserved_aas.second_cys, 96);
        assert_eq!(conserved_aas.j_trp_or_phe, 102);
    }
}
