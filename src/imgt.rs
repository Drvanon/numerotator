pub struct ConservedAminoAcids {
    first_cys: usize,
    conserved_trp: usize,
    hydrophobic_89: usize,
    second_cys: usize,
    j_trp_or_phe: usize,
}

use bio::alignment::{Alignment, AlignmentOperation};
use itertools::Itertools;

pub fn count_gaps_in_sequence_before_index(sequence: &[u8], index: usize) -> usize {
    sequence
        .into_iter()
        .take(index)
        .filter(|char| **char == b'-')
        .count()
}

#[derive(Debug)]
pub enum IMGTError {
    InvalidAlignment,
}

fn find_corresponding_position_in_alignment(alignment: &Alignment, pos: usize) -> Option<usize> {
    alignment
        .path()
        .into_iter()
        .find(|(x, y, op)| {
            *x == pos && (*op == AlignmentOperation::Match || *op == AlignmentOperation::Subst)
        })
        .map(|(_, y, _)| y)
}

#[derive(Debug)]
pub enum TransferErr {
    CouldNotFindPositionInAlignment,
}

impl ConservedAminoAcids {
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

        aa_23 == b'W'
            && aa_41 == b'C'
            && aa_104 == b'C'
            && [b'F', b'W'].contains(&aa_118)
            && [b'A', b'I', b'L', b'M', b'F', b'W', b'Y', b'V'].contains(&aa_89)
    }

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

    pub fn transfer(&self, alignment: &Alignment, destination: &[u8]) -> Result<Self, TransferErr> {
        Ok(Self {
            first_cys: find_corresponding_position_in_alignment(alignment, self.first_cys)
                .ok_or(TransferErr::CouldNotFindPositionInAlignment)?,
            conserved_trp: find_corresponding_position_in_alignment(alignment, self.conserved_trp)
                .ok_or(TransferErr::CouldNotFindPositionInAlignment)?,
            hydrophobic_89: find_corresponding_position_in_alignment(
                alignment,
                self.hydrophobic_89,
            )
            .ok_or(TransferErr::CouldNotFindPositionInAlignment)?,
            second_cys: find_corresponding_position_in_alignment(alignment, self.second_cys)
                .ok_or(TransferErr::CouldNotFindPositionInAlignment)?,
            j_trp_or_phe: find_corresponding_position_in_alignment(alignment, self.j_trp_or_phe)
                .ok_or(TransferErr::CouldNotFindPositionInAlignment)?,
        })
    }
}

#[cfg(test)]
mod test {
    use crate::initialize_ref_seqs;

    use super::*;

    #[test]
    fn test_validity_of_references() {
        let ref_seq = initialize_ref_seqs();
        ref_seq
            .into_iter()
            .for_each(|rec| assert!(ConservedAminoAcids::is_valid_alignment(rec.seq())))
    }

    #[test]
    fn test_imgtconserved_amino_acids_from_str() {
        let test_str = "QVQLVQSGA-EVKKPGASVKVSCKASGYTF----TSYGISWVRQAPGQGLEWMGWISAY--NGNTNYAQKLQ-GRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAR--------MDVWGQGTTVTVSS";
        let conserved_aas = ConservedAminoAcids::from_alignment(test_str.as_bytes()).unwrap();
        assert_eq!(conserved_aas.first_cys, 22);
        assert_eq!(conserved_aas.conserved_trp, 36);
        assert_eq!(conserved_aas.hydrophobic_89, 81);
        assert_eq!(conserved_aas.second_cys, 96);
        assert_eq!(conserved_aas.j_trp_or_phe, 102);
    }
}