#[derive(Clone, Debug)]
pub struct ReferenceSequence {
    alignment_str: String,
    name: String,
}

impl ReferenceSequence {
    pub fn get_sequence(&self) -> String {
        self.alignment_str.chars().filter(|c| *c != '-').collect()
    }

    pub fn get_missing_positions_in_framework(&self) {
        todo!()
    }
}
