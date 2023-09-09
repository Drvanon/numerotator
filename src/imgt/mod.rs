use std::collections::HashMap;

use thiserror::Error;
use tracing::trace;

use bio::{alignment::Alignment, io::fasta};

use self::{conserved_residues::ConservedResidues, reference::ReferenceSequence};

pub mod annotations;
pub mod conserved_residues;
pub mod reference;
pub mod regions;
pub mod single_letter_annotations;

const FR1_START: usize = 1;
const CDR1_START: usize = 27;
const FR2_START: usize = 39;
const CDR2_START: usize = 56;
const FR3_START: usize = 66;
const CDR3_START: usize = 105;
const FR4_START: usize = 118;
const FR4_END: usize = 128;

/// Error for when
#[derive(Debug, Error)]
pub enum IMGTError {
    #[error("Alignment did not have conserved residues in expected places.")]
    InvalidAlignment,

    #[error("Unexpected length ({1}) for region '{0}'.")]
    RegionTooLong(String, usize),

    #[error("CDR3 region too short. Expected at least 5, got {0}")]
    CDR3TooShort(usize),

    #[error("Region '{0}' and '{0}' overlapped.")]
    OverlappingRegions(String, String),

    #[error("Bad alignment string.")]
    BadBytesInAlignment(#[from] std::str::Utf8Error),
}

/// Error thrown when looking for a reference sequence.
#[derive(Debug, Error)]
pub enum RefSeqErr {
    #[error("Could not find reference record for record {0}")]
    NoReferenceSequenceFound(fasta::Record),
}

/// Captures an alignment of a query sequence to reference sequence.
///
/// Uses records to keep track of identities. (For the reference this
/// is particularly important since often you will want to look up the
/// associated curated alignment.)
pub struct ReferenceAlignment {
    pub reference: ReferenceSequence,
    pub query_record: fasta::Record,
    pub alignment: Alignment,
}

/// Find the record that produces the best alignment.
pub fn find_best_reference_sequence(
    record: fasta::Record,
    ref_seqs: &HashMap<&str, ReferenceSequence>,
) -> Result<ReferenceAlignment, RefSeqErr> {
    trace!(query_seq = record.id(), "Finding reference sequence.");
    // TODO: Optimize settings.
    // Settings taken from rust bio example. Fully unoptimized.
    let mut aligner =
        bio::alignment::pairwise::Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });

    // TODO: Optimize this to go by alignment block!
    ref_seqs
        .values()
        .map(|reference_sequence| {
            (
                reference_sequence,
                aligner.local(&reference_sequence.get_sequence(), record.seq()),
            )
        })
        .max_by_key(|(_reference, alignment)| alignment.score)
        .map(|(reference, alignment)| {
            trace!(
                score = alignment.score,
                reference = reference.name,
                "Found alignment."
            );
            ReferenceAlignment {
                // Cloning here should not be a huge problem, since we only clone once per query sequence.
                reference: reference.clone(),
                alignment,
                query_record: record.clone(),
            }
        })
        .ok_or(RefSeqErr::NoReferenceSequenceFound(record))
}
