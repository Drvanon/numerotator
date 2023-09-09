use bio::alignment::AlignmentOperation;

/// Numbering of single amino acids.
///
/// Mapping according to [this](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html) IMGT scientific chart.
use super::annotations::{Annotation, VRegionAnnotation};
use super::{IMGTError, ReferenceAlignment};
use crate::imgt;
use std::collections::HashMap;

fn number_cdr1(start: usize, end: usize) -> Result<Vec<Annotation>, IMGTError> {
    let cdr1_size = end - start;
    let cdr1_length_ranges_mapping: HashMap<usize, Vec<usize>> = [
        (12, vec![27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38]),
        (11, vec![27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38]),
        (10, vec![27, 28, 29, 30, 31, 34, 35, 36, 37, 38]),
        (9, vec![27, 28, 29, 30, 31, 35, 36, 37, 38]),
        (8, vec![27, 28, 29, 30, 35, 36, 37, 38]),
        (7, vec![27, 28, 29, 30, 36, 37, 38]),
        (6, vec![27, 28, 29, 36, 37, 38]),
        (5, vec![27, 28, 29, 37, 38]),
    ]
    .into_iter()
    .collect();

    Ok(cdr1_length_ranges_mapping
        .get(&cdr1_size)
        .ok_or(IMGTError::RegionTooLong("CDR1-IMGT".to_string(), cdr1_size))?
        .into_iter()
        .map(|number| number.to_string())
        .zip(start..end)
        .map(|(name, position)| Annotation {
            start: position,
            end: position + 1,
            name,
        })
        .collect())
}

fn number_cdr2(start: usize, end: usize) -> Result<Vec<Annotation>, IMGTError> {
    let cdr2_size = end - start;
    let cdr2_length_ranges_mapping: HashMap<usize, Vec<usize>> = [
        (10, vec![56, 57, 58, 59, 60, 61, 62, 63, 64, 65]),
        (9, vec![56, 57, 58, 59, 60, 62, 63, 64, 65]),
        (8, vec![56, 57, 58, 59, 62, 63, 64, 65]),
        (7, vec![56, 57, 58, 59, 63, 64, 65]),
        (6, vec![56, 57, 58, 63, 64, 65]),
        (5, vec![56, 57, 58, 64, 65]),
        (4, vec![56, 57, 64, 65]),
        (3, vec![56, 57, 65]),
        (2, vec![56, 65]),
        (1, vec![56]),
        (0, vec![]),
    ]
    .into_iter()
    .collect();

    Ok(cdr2_length_ranges_mapping
        .get(&cdr2_size)
        .ok_or(IMGTError::RegionTooLong("CDR2-IMGT".to_string(), cdr2_size))?
        .into_iter()
        .map(|number| number.to_string())
        .zip(start..end)
        .map(|(name, position)| Annotation {
            start: position,
            end: position + 1,
            name,
        })
        .collect())
}

fn number_cdr3(start: usize, end: usize) -> Result<Vec<Annotation>, IMGTError> {
    let cdr3_size = end - start;
    if cdr3_size < 5 {
        return Err(IMGTError::CDR3TooShort(cdr3_size));
    }

    let cdr3_length_ranges_mapping: HashMap<usize, Vec<usize>> = [
        (
            13,
            vec![
                105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
            ],
        ),
        (
            12,
            vec![105, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117],
        ),
        (
            11,
            vec![105, 106, 107, 108, 109, 110, 113, 114, 115, 116, 117],
        ),
        (10, vec![105, 106, 107, 108, 109, 113, 114, 115, 116, 117]),
        (9, vec![105, 106, 107, 108, 109, 114, 115, 116, 117]),
        (8, vec![105, 106, 107, 108, 114, 115, 116, 117]),
        (7, vec![105, 106, 107, 108, 115, 116, 117]),
        (6, vec![105, 106, 107, 115, 116, 117]),
        (5, vec![105, 106, 107, 116, 117]),
    ]
    .into_iter()
    .collect();

    if cdr3_size <= 13 {
        return Ok(cdr3_length_ranges_mapping
            .get(&cdr3_size)
            .unwrap()
            .into_iter()
            .map(|number| number.to_string())
            .zip(start..end)
            .map(|(name, position)| Annotation {
                start: position,
                end: position + 1,
                name,
            })
            .collect());
    }

    // TODO: extract this range-zip-map stuff to a function for DRY.
    Ok((imgt::CDR3_START..=111)
        .zip(start..start + 6)
        .map(|(number, position)| Annotation {
            start: position,
            end: position + 1,
            name: number.to_string(),
        })
        .chain(additional_positions_between_111_and_112(start + 6, end - 5).into_iter())
        .chain(
            (113..imgt::FR4_START)
                .zip(end - 5..end)
                .map(|(number, position)| Annotation {
                    start: position,
                    end: position + 1,
                    name: number.to_string(),
                }),
        )
        .collect())
}

/// Additional positions between 111 and 112 in the CDR3-IMGT region.
fn additional_positions_between_111_and_112(start: usize, end: usize) -> Vec<Annotation> {
    let n_extra_positions = end - start;
    let n_extra_positions_111 = (n_extra_positions as f64 / 2.0).floor() as usize;
    let n_extra_positions_112 = (n_extra_positions as f64 / 2.0).ceil() as usize;

    let extra_positions_111 = (0..n_extra_positions_111)
        .into_iter()
        .map(|i| format!("111.{}", i));
    let extra_positions_112 = (0..n_extra_positions_112)
        .into_iter()
        .map(|i| format!("112.{}", i))
        .rev();

    extra_positions_111
        .chain(extra_positions_112)
        .zip(start..end)
        .map(|(name, position)| Annotation {
            start: position,
            end: position + 1,
            name,
        })
        .collect()
}

fn number_framework(
    reference_alignment: &ReferenceAlignment,
    framework: imgt::Framework,
) -> Vec<Annotation> {
    let range = match framework {
        imgt::Framework::FR1 => imgt::FR1,
        imgt::Framework::FR2 => imgt::FR2,
        imgt::Framework::FR3 => imgt::FR3,
        imgt::Framework::FR4 => imgt::FR4,
    };
    let numbers = range.clone().filter(|pos| {
        reference_alignment
            .reference
            .get_missing_positions_in_framework(&framework)
            .contains(pos)
    });
    reference_alignment
        .alignment
        .path()
        .into_iter()
        .filter(|(x, _y, _op)| range.contains(&x))
        .flat_map(|(_x, y, op)| match op {
            AlignmentOperation::Match => Some(y),
            AlignmentOperation::Subst => Some(y),
            AlignmentOperation::Del => None,
            AlignmentOperation::Ins => None,
            AlignmentOperation::Xclip(_) => None,
            AlignmentOperation::Yclip(_) => None,
        })
        .zip(numbers)
        .map(|(number, position)| Annotation {
            // Path starts at one, where as annotations are zero based.
            start: position - 1,
            end: position,
            name: number.to_string(),
        })
        .collect()
}

impl VRegionAnnotation {
    pub fn number_regions(
        &self,
        reference_alignment: &ReferenceAlignment,
    ) -> Result<Vec<Annotation>, IMGTError> {
        Ok(number_framework(&reference_alignment, imgt::Framework::FR1)
            .into_iter()
            .chain(
                number_cdr1(self.cdr_annotation.cdr1.start, self.cdr_annotation.cdr1.end)?
                    .into_iter(),
            )
            .chain(number_framework(&reference_alignment, imgt::Framework::FR2))
            .chain(
                number_cdr2(self.cdr_annotation.cdr2.start, self.cdr_annotation.cdr2.end)?
                    .into_iter(),
            )
            .chain(number_framework(&reference_alignment, imgt::Framework::FR3))
            .chain(
                number_cdr3(self.cdr_annotation.cdr3.start, self.cdr_annotation.cdr3.end)?
                    .into_iter(),
            )
            .chain(number_framework(&reference_alignment, imgt::Framework::FR4))
            .collect())
    }
}
