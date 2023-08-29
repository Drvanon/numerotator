use bio::{
    alignment::{self, Alignment},
    io::fasta,
};
use clap::{value_parser, Parser};
use numerator::imgt;
use std::path::PathBuf;
use thiserror::Error;
use tracing::{debug, error, info, trace, Level};
use tracing_subscriber::FmtSubscriber;

#[derive(Parser, Debug)]
#[command()]
struct Args {
    #[arg(index = 1, num_args=..)]
    sequences: Vec<String>,
    #[arg(short, long, value_parser=value_parser!(PathBuf))]
    sequences_file: Option<PathBuf>,

    #[arg(short, long)]
    only_regions: bool,
}

#[derive(Debug, Error)]
pub enum RefSeqErr {
    #[error("Could not find reference record for record {0}")]
    NoReferenceSequenceFound(fasta::Record),
}

struct ReferenceAlignment {
    reference_record: fasta::Record,
    query_record: fasta::Record,
    alignment: Alignment,
}

fn find_best_reference_sequence(
    record: fasta::Record,
    ref_seqs: &Vec<fasta::Record>,
) -> Result<ReferenceAlignment, RefSeqErr> {
    let mut aligner =
        alignment::pairwise::Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });
    // TODO: Optimize this to go by alignment block!
    ref_seqs
        .into_iter()
        .map(|reference_record| {
            (
                reference_record,
                aligner.local(record.seq(), reference_record.seq()),
            )
        })
        .max_by_key(|(_reference, alignment)| alignment.score)
        .map(|(reference_record, alignment)| {
            trace!(
                score = alignment.score,
                reference = reference_record.id(),
                "Found alignment."
            );
            ReferenceAlignment {
                // Cloning here should not be a huge problem, since we only clone once per query sequence.
                reference_record: reference_record.clone(),
                alignment,
                query_record: record.clone(),
            }
        })
        .ok_or(RefSeqErr::NoReferenceSequenceFound(record))
}

fn report_error<OkType, ErrType: std::fmt::Display>(
    result: Result<OkType, ErrType>,
) -> Result<OkType, ErrType> {
    result.map_err(|err| {
        error!("{}", err);
        err
    })
}
fn transfer_conserved_residues_via_alignment(
    reference_alignment: ReferenceAlignment,
    all_conserved_residues: &std::collections::HashMap<&str, imgt::ConservedResidues>,
) -> Result<(imgt::ConservedResidues, ReferenceAlignment), imgt::TransferErr> {
    trace!(
        reference = reference_alignment.reference_record.id(),
        "Looking for reference alignment"
    );
    all_conserved_residues
        .get(reference_alignment.reference_record.id())
        .expect("Reference sequence id should be in reference alignments aswell.")
        .transfer(
            &reference_alignment.alignment,
            reference_alignment.query_record.seq(),
        )
        .and_then(|conserved_residues| Ok((conserved_residues, reference_alignment)))
}

fn main() {
    let args = Args::parse();

    let subscriber = FmtSubscriber::builder()
        // all spans/events with a level higher than TRACE (e.g, debug, info, warn, etc.)
        .with_max_level(Level::TRACE)
        // will be written to stderr
        .with_writer(std::io::stderr)
        // completes the builder.
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("setting default subscriber failed");

    debug!("Initializing reference sequences.");
    let ref_seqs = imgt::initialize_ref_seqs();

    debug!("Initializing conserved residues.");
    let all_conserved_residues = imgt::initialize_conserved_residues();

    debug!("Collecting sequences from command line.");
    let sequences_from_command_line = args.sequences.into_iter().enumerate().map(|(i, seq)| {
        fasta::Record::with_attrs(
            i.to_string().as_str(), // TODO: Use uuid here in order to prevent clash with potential use case.
            Some(format!("sequence {} from the command line", i).as_str()),
            seq[..].as_bytes(),
        )
    });

    let sequences_from_sequence_file = args.sequences_file.and_then(|path| {
        info!("Reading input sequences file.");
        Some(
            fasta::Reader::new(std::fs::File::open(path).expect("Could not open sequences file."))
                .records()
                .map(|record_result| {
                    record_result.expect("Could not parse record in sequences file.")
                }),
        )
    });

    let vregion_annotations: Vec<(_, _)> = sequences_from_command_line
        .chain(sequences_from_sequence_file.into_iter().flatten())
        .map(|query_seq| {
            trace!(query_seq = query_seq.id(), "Finding reference sequence.");
            find_best_reference_sequence(query_seq, &ref_seqs)
        })
        .flat_map(report_error)
        .map(|reference_alignment| {
            transfer_conserved_residues_via_alignment(reference_alignment, &all_conserved_residues)
        })
        .flat_map(report_error)
        .map(
            |(conserved_residues, reference_alignment)| -> Result<_, imgt::annotations::AnnotationError> {
                trace!(query_seq = reference_alignment.query_record.id(), "Creating VREGION annotation.");
                Ok((
                    imgt::annotations::VRegionAnnotation::try_from(conserved_residues)?,
                    reference_alignment,
                ))
            },
        )
        .flat_map(report_error)
        .collect();

    if args.only_regions {
        let mut writer = fasta::Writer::new(std::io::stdout());

        vregion_annotations
            .into_iter()
            .map(|(vregion_annotation, reference_alignment)| -> Vec<_> {
                vregion_annotation
                    .into_iter()
                    .map(|ann| {
                        imgt::annotations::apply_annotation(&reference_alignment.query_record, &ann)
                    })
                    .collect()
            })
            .flatten()
            .for_each(|record| {
                writer
                    .write_record(&record)
                    .expect("Could not write record.")
            });
    }
}
