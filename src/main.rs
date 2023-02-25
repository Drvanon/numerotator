use std::{collections::HashMap, ops::RangeInclusive};

use itertools::Itertools;

static NC10_SEQ: &str = "QVQLQQSGAELVKPGASVRMSCKASGYTFTNYNMYWVKQSPGQGLEWIGIFYPGNGDTSYNQKFKDKATLTADKSSNTAYMQLSSLTSEDSAVYYCARSGGSYRYDGGFDYWGQGTTVTV";
static TEST_SEQ: &str = "AAAAAAQVQLQQSGAELVKPGASVRMSCKASGYTFTNYNAAAAAAAAMYWVKQSPGQGLEWIGIFYPGNGDTSYNQKFKDKATLTADKSSNTAYMQLSSLTSEAAAAAAAAADSAVYYCARSGGSYRYDGGFDYWGQGTTVTV";

type Kmer = String;
type Name = String;
type Position = usize;

struct SequenceProfile(HashMap<Kmer, Position>);
impl From<&str> for SequenceProfile {
    fn from(value: &str) -> Self {
        let v: Vec<char> = value.chars().collect();
        Self(
            v.windows(7)
                .enumerate()
                .map(|(i, w)| (w.iter().collect(), i))
                .collect(),
        )
    }
}

type DiagonalIndex = isize;

#[derive(Debug, Default)]
struct KmerHits(Vec<(Position, Position)>);
impl KmerHits {
    fn find_ranges(numbers: Vec<usize>) -> Vec<RangeInclusive<usize>> {
        let gaps: Vec<(usize, usize)> = numbers
            .iter()
            .sorted()
            .tuple_windows()
            .filter(|(prev, next)| (*prev + 1) == **next)
            .map(|(p, n)| (*p, *n))
            .collect();

        let mut ranges: Vec<RangeInclusive<usize>> = gaps
            .iter()
            .tuple_windows()
            .map(|(prev_gap, next_gap)| RangeInclusive::new(prev_gap.1, next_gap.0))
            .collect();

        if let Some(last_gap) = gaps.last() {
            if let Some(&last_number) = numbers.last() {
                ranges.push(RangeInclusive::new(last_gap.1, last_number));
            }
        }

        if let Some(last_gap) = gaps.first() {
            if let Some(&last_number) = numbers.last() {
                ranges.push(RangeInclusive::new(last_gap.1, last_number));
            }
        }

        ranges
    }

    pub fn get_identical_kmer_ranges(&self) -> Vec<(DiagonalIndex, Vec<RangeInclusive<usize>>)> {
        let diagonals: HashMap<DiagonalIndex, Vec<Position>> = self
            .0
            .iter()
            .map(|(pos, ref_pos)| (*pos as isize - *ref_pos as isize, *pos))
            .into_group_map();
        diagonals
            .into_iter()
            .map(|(diagonal_index, kmer_starts)| (diagonal_index, Self::find_ranges(kmer_starts)))
            .collect()
    }
}

#[derive(Debug)]
struct KmerDatabase(HashMap<Kmer, Vec<(Name, Position)>>);
impl Default for KmerDatabase {
    fn default() -> Self {
        Self(Default::default())
    }
}
impl KmerDatabase {
    pub fn add_profile(&mut self, name: String, profile: SequenceProfile) {
        profile.0.into_iter().for_each(|(kmer, position)| {
            self.0
                .entry(kmer)
                .or_insert(vec![])
                .push((name.clone(), position))
        });
    }

    pub fn find_kmer_hits(&self, profile: SequenceProfile) -> HashMap<Name, KmerHits> {
        profile
            .0
            .into_iter()
            .flat_map(|(kmer, pos)| {
                self.0.get(&kmer).map(move |results| {
                    results
                        .into_iter()
                        .map(move |(name, ref_position)| (name.clone(), (pos, *ref_position)))
                })
            })
            .flatten()
            .into_group_map()
            .into_iter()
            .map(|(name, hits)| (name, KmerHits(hits)))
            .collect()
    }
}

fn main() {
    let mut database = KmerDatabase::default();
    database.add_profile("NC10".into(), NC10_SEQ.into());
    println!("{:?}", database.find_kmer_hits(TEST_SEQ.into()));
}
