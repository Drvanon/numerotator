use bio::{
    alignment::{self, Alignment},
    io::fasta,
};
use clap::{value_parser, Parser};
use numerator::imgt::{self, annotations::VRegionAnnotation};
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
struct ReferenceAlignment {
    reference_record: fasta::Record,
    query_record: fasta::Record,
    alignment: Alignment,
}

/// Find the record that produces the best alignment.
fn find_best_reference_sequence(
    record: fasta::Record,
    ref_seqs: &Vec<fasta::Record>,
) -> Result<ReferenceAlignment, RefSeqErr> {
    // TODO: Optimize settings.
    // Settings taken from rust bio example. Fully unoptimized.
    let mut aligner =
        alignment::pairwise::Aligner::new(-5, -1, |a, b| if a == b { 1i32 } else { -1i32 });

    // TODO: Optimize this to go by alignment block!
    ref_seqs
        .into_iter()
        .map(|reference_record| {
            (
                reference_record,
                aligner.local(reference_record.seq(), record.seq()),
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

    info!("Initializing...");
    debug!("Initializing reference sequences.");
    let ref_seqs = imgt::initialize_ref_seqs();

    debug!("Initializing conserved residues.");
    let all_conserved_residues = imgt::initialize_conserved_residues();

    // Records are much nicer to deal with than straigt strings, since they carry their own
    // identifier and description. Now they don't have to be generated at the call site.
    // It might not be great to be tied to fasta though.
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

    sequences_from_command_line
        .chain(sequences_from_sequence_file.into_iter().flatten())
        .map(|query_seq| {
            trace!(query_seq = query_seq.id(), "Finding reference sequence.");
            find_best_reference_sequence(query_seq, &ref_seqs)
        })
        .flat_map(report_error)
        .map(|reference_alignment| -> Result<(VRegionAnnotation, ReferenceAlignment), anyhow::Error> {
            let reference_conserved_residues = all_conserved_residues
                .get(reference_alignment.reference_record.id())
                .expect("Reference sequence id should be in reference alignments aswell.");
            trace!(
                query_seq = reference_alignment.query_record.id(),
                alignment = format!("{:?}", reference_alignment.alignment.path()),
                "Transferring reference alignment."
            );
            let vregions = reference_conserved_residues
                .transfer(
                    &reference_alignment.alignment,
                    reference_alignment.query_record.seq(),
                )
                .map_err(anyhow::Error::from)
                .and_then(|conserved_residues| {
                    trace!(
                        query_seq = reference_alignment.query_record.id(),
                        "Creating VREGION annotation."
                    );
                    imgt::annotations::VRegionAnnotation::try_from(
                        &conserved_residues,
                        &reference_alignment.alignment,
                    )
                    .map_err(anyhow::Error::from)
                });
            Ok((vregions?, reference_alignment))
        })
        .flat_map(report_error)
        .for_each(|(vregion_annotation, reference_alignment)| {
            if args.only_regions {
                trace!(
                    query_seq = reference_alignment.query_record.id(),
                    "Applying annotations."
                );
                write_vregion_annotations(
                    &reference_alignment.query_record,
                    &vregion_annotation,
                    std::io::stderr(),
                )
            }
        });
}

/// Apply all annotations of the a vregion to a record and write them to a writer.
fn write_vregion_annotations<W: std::io::Write>(
    record: &fasta::Record,
    vregion_annotation: &imgt::annotations::VRegionAnnotation,
    writer: W,
) {
    let mut fasta_writer = fasta::Writer::new(writer);
    vregion_annotation
        .clone()
        .into_iter()
        .map(|ann| imgt::annotations::apply_annotation(&record, &ann))
        .for_each(|record| {
            fasta_writer
                .write_record(&record)
                .expect("Could not write record.")
        });
}
