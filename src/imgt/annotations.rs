use bio::io::fasta;
use thiserror::Error;

use super::ConservedResidues;

#[derive(Clone)]
pub struct Annotation {
    pub start: usize,
    pub end: usize,
    pub name: String,
}

pub fn apply_annotation(record: &fasta::Record, annotation: &Annotation) -> fasta::Record {
    fasta::Record::with_attrs(
        format!("{}_{}", record.id(), annotation.name).as_str(),
        Some(format!("{} of {}", annotation.name, record.id()).as_str()),
        &record.seq()[annotation.start..annotation.end],
    )
}

#[derive(Clone)]
pub struct FrameworkAnnotation {
    pub fr1: Annotation,
    pub fr2: Annotation,
    pub fr3: Annotation,
    pub fr4: Annotation,
}

impl IntoIterator for FrameworkAnnotation {
    type Item = Annotation;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        vec![self.fr1, self.fr2, self.fr3, self.fr4].into_iter()
    }
}

#[derive(Clone)]
pub struct CDRAnnotation {
    pub cdr1: Annotation,
    pub cdr2: Annotation,
    pub cdr3: Annotation,
}
impl IntoIterator for CDRAnnotation {
    type Item = Annotation;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        vec![self.cdr1, self.cdr2, self.cdr3].into_iter()
    }
}

#[derive(Clone)]
pub struct VRegionAnnotation {
    pub cdr_annotation: CDRAnnotation,
    pub framework_annotation: FrameworkAnnotation,
}

impl IntoIterator for VRegionAnnotation {
    type Item = Annotation;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        vec![
            self.framework_annotation.fr1,
            self.cdr_annotation.cdr1,
            self.framework_annotation.fr2,
            self.cdr_annotation.cdr2,
            self.framework_annotation.fr3,
            self.cdr_annotation.cdr3,
            self.framework_annotation.fr4,
        ]
        .into_iter()
    }
}

#[derive(Debug, Error, Clone)]
pub enum AnnotationError {
    #[error("Region '{0}' and '{0}' overlapped.")]
    OverlappingRegions(String, String),
}

impl TryFrom<ConservedResidues> for FrameworkAnnotation {
    type Error = AnnotationError;

    fn try_from(conserved_residues: ConservedResidues) -> Result<Self, AnnotationError> {
        // TODO: this does not take into account a case where there is a G missing in the GGG sequence.
        let fr1 = Annotation {
            start: conserved_residues.first_cys.saturating_sub(23),
            end: conserved_residues.first_cys + 3,
            name: "FR1-IMGT".to_string(),
        };
        let fr2 = Annotation {
            start: conserved_residues.conserved_trp - 2,
            end: conserved_residues.conserved_trp + 14,
            name: "FR2-IMGT".to_string(),
        };
        let fr3 = Annotation {
            start: conserved_residues.hydrophobic_89 - 23,
            end: conserved_residues.second_cys,
            name: "FR3-IMGT".to_string(),
        };
        let fr4 = Annotation {
            start: conserved_residues.j_trp_or_phe,
            end: conserved_residues.j_trp_or_phe + 10,
            name: "FR4-IMGT".to_string(),
        };

        if fr1.end > fr2.start {
            return Err(AnnotationError::OverlappingRegions(fr1.name, fr2.name));
        }

        if fr2.end > fr3.start {
            return Err(AnnotationError::OverlappingRegions(fr2.name, fr3.name));
        }

        if fr3.end > fr4.start {
            return Err(AnnotationError::OverlappingRegions(fr3.name, fr4.name));
        }

        Ok(Self { fr1, fr2, fr3, fr4 })
    }
}

impl TryFrom<FrameworkAnnotation> for CDRAnnotation {
    type Error = AnnotationError;
    fn try_from(framework_annotation: FrameworkAnnotation) -> Result<Self, Self::Error> {
        let cdr1 = Annotation {
            start: framework_annotation.fr1.end,
            end: framework_annotation.fr2.start,
            name: "CDR1-IMGT".to_string(),
        };

        let cdr2 = Annotation {
            start: framework_annotation.fr2.end,
            end: framework_annotation.fr3.start,
            name: "CDR2-IMGT".to_string(),
        };

        let cdr3 = Annotation {
            start: framework_annotation.fr3.end,
            end: framework_annotation.fr4.start,
            name: "CDR3-IMGT".to_string(),
        };

        Ok(Self { cdr1, cdr2, cdr3 })
    }
}

impl TryFrom<ConservedResidues> for VRegionAnnotation {
    type Error = AnnotationError;

    fn try_from(conserved_residues: ConservedResidues) -> Result<Self, Self::Error> {
        let framework_annotation = FrameworkAnnotation::try_from(conserved_residues)?;
        let cdr_annotation = CDRAnnotation::try_from(framework_annotation.clone())?;
        Ok(Self {
            framework_annotation,
            cdr_annotation,
        })
    }
}
