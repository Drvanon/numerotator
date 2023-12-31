use bio::io::fasta;

/// Annotation of a sequence.
#[derive(Clone)]
pub struct Annotation {
    pub start: usize,
    pub end: usize,
    pub name: String,
}

/// Create a new record for the subsequence that the annotation references in a given record.
pub fn apply_annotation(record: &fasta::Record, annotation: &Annotation) -> fasta::Record {
    fasta::Record::with_attrs(
        format!("{}_{}", annotation.name, record.id()).as_str(),
        Some(
            format!(
                "IMGT Number {} on {}|{}|{}",
                annotation.name,
                record.id(),
                annotation.start,
                annotation.end
            )
            .as_str(),
        ),
        &record.seq()[annotation.start..annotation.end],
    )
}

/// IMGT Framework (FRx-IMGT) annotations of a VREGION sequence.
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

/// IMGT CDR (CDRx-IMGT) annotations of a VREGION sequence.
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

/// VREGION annotation of a sequence.
#[derive(Clone)]
pub struct VRegionAnnotation {
    pub cdr_annotation: CDRAnnotation,
    pub framework_annotation: FrameworkAnnotation,
}

impl VRegionAnnotation {
    pub fn region_annotations(&self) -> Vec<Annotation> {
        vec![
            self.framework_annotation.fr1.clone(),
            self.cdr_annotation.cdr1.clone(),
            self.framework_annotation.fr2.clone(),
            self.cdr_annotation.cdr2.clone(),
            self.framework_annotation.fr3.clone(),
            self.cdr_annotation.cdr3.clone(),
            self.framework_annotation.fr4.clone(),
        ]
    }
}
