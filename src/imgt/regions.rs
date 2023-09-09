use bio::alignment::{Alignment, AlignmentOperation};

use super::{
    annotations::{Annotation, CDRAnnotation, FrameworkAnnotation, VRegionAnnotation},
    ConservedResidues, IMGTError,
};

// TODO: Find a better name than try from.
// TODO: Should there be an option to not use an alignment?

impl FrameworkAnnotation {
    /// Try to create the framework annotations.
    ///
    /// Alignment assumes that sequence x was an IMGT reference sequence.
    fn try_from(
        conserved_residues: &ConservedResidues,
        alignment: &Alignment,
    ) -> Result<Self, IMGTError> {
        let v_region_start_position = alignment
            .path()
            .into_iter()
            .find(|(x, _, op)| {
                *x == 1
                    && match op {
                        AlignmentOperation::Xclip(_) => false,
                        AlignmentOperation::Yclip(_) => false,
                        _ => true,
                    }
            })
            .expect("Alignment should contain first amino acid of reference sequence.")
            .1;

        // bio::alignment::Alignment::path uses 1 based indexing.
        let v_region_start = v_region_start_position - 1;

        let v_region_end = alignment
            .path()
            .into_iter()
            .find(|(x, _, op)| {
                *x == alignment.xend
                    && match op {
                        AlignmentOperation::Xclip(_) => false,
                        AlignmentOperation::Yclip(_) => false,
                        _ => true,
                    }
            })
            .expect("Alignment should contain first amino acid of reference sequence.")
            .1;

        let fr1 = Annotation {
            start: v_region_start,
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
            end: v_region_end,
            name: "FR4-IMGT".to_string(),
        };

        if fr1.end > fr2.start {
            return Err(IMGTError::OverlappingRegions(fr1.name, fr2.name));
        }

        if fr2.end > fr3.start {
            return Err(IMGTError::OverlappingRegions(fr2.name, fr3.name));
        };

        if fr3.end > fr4.start {
            return Err(IMGTError::OverlappingRegions(fr3.name, fr4.name));
        }

        Ok(Self { fr1, fr2, fr3, fr4 })
    }

    pub fn get_cdr1_length(&self) -> usize {
        self.fr2.start - self.fr1.end
    }
    pub fn get_cdr2_length(&self) -> usize {
        self.fr3.start - self.fr2.end
    }
    pub fn get_cdr3_length(&self) -> usize {
        self.fr4.start - self.fr3.end
    }
}

impl TryFrom<FrameworkAnnotation> for CDRAnnotation {
    type Error = IMGTError;

    /// Try to create the framework annotations.
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

impl VRegionAnnotation {
    /// Try to create a VREGION annotation from the positions of conserved residues and an alignment.
    pub fn try_from(
        conserved_residues: &ConservedResidues,
        alignment: &Alignment,
    ) -> Result<Self, IMGTError> {
        let framework_annotation = FrameworkAnnotation::try_from(conserved_residues, alignment)?;
        let cdr_annotation = CDRAnnotation::try_from(framework_annotation.clone())?;
        Ok(Self {
            framework_annotation,
            cdr_annotation,
        })
    }
}
