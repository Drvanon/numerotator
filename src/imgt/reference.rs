use std::collections::HashMap;

use itertools::Itertools;

use super::{conserved_residues::ConservedResidues, IMGTError};

pub fn is_valid_alignment(alignment: &[u8]) -> Option<ConservedResidues> {
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
        None => return None,
    };

    if aa_23 == b'C'
        && aa_41 == b'W'
        && aa_104 == b'C'
        && [b'F', b'W'].contains(&aa_118)
        && [b'A', b'I', b'L', b'M', b'F', b'W', b'Y', b'V'].contains(&aa_89)
    {
        Some(ConservedResidues::from(alignment))
    } else {
        None
    }
}

#[derive(Clone, Debug)]
pub struct ReferenceSequence {
    alignment: String,
    pub name: String,
    conserved_residues: ConservedResidues,
}

impl ReferenceSequence {
    pub fn new(name: &str, alignment: &[u8]) -> Result<Self, IMGTError> {
        Ok(Self {
            alignment: std::str::from_utf8(alignment)
                .map_err(|e| IMGTError::from(e))?
                .to_string(),
            name: name.to_string(),
            conserved_residues: is_valid_alignment(alignment).ok_or(IMGTError::InvalidAlignment)?,
        })
    }

    pub fn get_conserved_residues(&self) -> &ConservedResidues {
        &self.conserved_residues
    }

    pub fn get_sequence(&self) -> Vec<u8> {
        self.alignment
            .as_bytes()
            .into_iter()
            .map(|&b| b)
            .filter(|c| *c != b'-')
            .collect()
    }

    pub fn get_missing_positions_in_fr1(&self) -> Vec<usize> {
        todo!()
    }

    pub fn get_alignment(&self) -> &[u8] {
        self.alignment.as_bytes()
    }
}

/// Load the precomputed and curated reference sequences.
pub fn initialize_reference_sequences() -> HashMap<&'static str, ReferenceSequence> {
    // TODO: Write a proper stockholm reader.
    let stockholm_data = include_str!("../data/reference.stockholm");
    stockholm_data
        .split_ascii_whitespace()
        .tuples()
        .filter_map(|(id, alignment)| {
            Some((id, ReferenceSequence::new(id, alignment.as_bytes()).ok()?))
        })
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;
    use tracing_test::traced_test;
    const TEST_ALIGNMENT_STR: &str = "QVQLVQSGA-EVKKPGASVKVSCKASGYTF----TSYGISWVRQAPGQGLEWMGWISAY--NGNTNYAQKLQ-GRVTMTTDTSTSTAYMELRSLRSDDTAVYYCAR--------MDVWGQGTTVTVSS";

    #[test]
    #[traced_test]
    fn test_validity() {
        assert!(is_valid_alignment(TEST_ALIGNMENT_STR.as_bytes()).is_some())
    }

    #[test]
    fn test_validity_of_references() {
        let ref_seqs = initialize_reference_sequences();
        ref_seqs
            .values()
            .for_each(|rec| assert!(is_valid_alignment(rec.get_alignment()).is_some()))
    }

    #[test]
    fn test_new_reference_sequence() {
        let ref_seq_res = ReferenceSequence::new("test", TEST_ALIGNMENT_STR.as_bytes());
        assert!(ref_seq_res.is_ok());
        assert_eq!(ref_seq_res.unwrap().name, "test");
    }
}
