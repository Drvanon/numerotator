use bio::io::fasta;
use clap::{value_parser, Parser};
use numerotator::imgt::{
    self,
    annotations::{Annotation, VRegionAnnotation},
    find_best_reference_sequence, ReferenceAlignment, conserved_residues::ConservedResidues,
};
use std::path::PathBuf;
use tracing::{debug, error, info, trace, Level};
use tracing_subscriber::FmtSubscriber;

#[derive(Parser, Debug)]
#[command()]
struct Args {
    #[arg(index = 1, num_args=..)]
    sequences: Vec<String>,
    #[arg(short, long, value_parser=value_parser!(PathBuf))]
    sequences_file: Option<PathBuf>,

    #[arg(short, long, help = "Annotate the regions as well.")]
    annotate_regions: bool,

    #[arg(
        short,
        long,
        help = "Do not number the sequences. (Useful in combination with --annotate-regions)"
    )]
    no_number: bool,
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
    let ref_seqs = imgt::reference::initialize_reference_sequences();

    // Records are much nicer to deal with than simple strings, since they carry their own
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
        .map(|query_seq| find_best_reference_sequence(query_seq, &ref_seqs) )
        .flat_map(report_error)
        .map(|reference_alignment| -> Result<(VRegionAnnotation, ReferenceAlignment), anyhow::Error> {
            trace!(
                query_seq = reference_alignment.query_record.id(),
                alignment = format!("{:?}", reference_alignment.alignment.path()),
                "Transferring reference alignment."
            );
            let vregions = transfer_conserved_residues(reference_alignment.reference.get_conserved_residues(), &reference_alignment);
            Ok((vregions?, reference_alignment))
        })
        .flat_map(report_error)
        .for_each(|(vregion_annotation, reference_alignment)| {
            if args.annotate_regions {
                trace!(
                    query_seq = reference_alignment.query_record.id(),
                    "Applying region annotations."
                );
                write_annotations(
                    &reference_alignment.query_record,
                    vregion_annotation.region_annotations(),
                    std::io::stdout(),
                );
            }

            if !args.no_number {
                trace!("Applying numbering.");
                let number_annotations =  vregion_annotation.number_regions(&reference_alignment);
                match number_annotations {
                    Ok(annotations) => {
                        write_annotations(&reference_alignment.query_record, annotations , std::io::stdout())
                    },
                    Err(error) => {
                        error!(sequence = reference_alignment.query_record.id(), error=error.to_string(), "Could not number regions for sequence.");
                    }
                }
                                
            }
        });
}

fn transfer_conserved_residues(
    reference_conserved_residues: &ConservedResidues,
    reference_alignment: &ReferenceAlignment,
) -> Result<VRegionAnnotation, anyhow::Error> {
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
    vregions
}

/// Apply all annotations of the a vregion to a record and write them to a writer.
fn write_annotations<W: std::io::Write>(
    record: &fasta::Record,
    annotations: Vec<Annotation>,
    writer: W,
) {
    let mut fasta_writer = fasta::Writer::new(writer);
    annotations
        .into_iter()
        .map(|ann| imgt::annotations::apply_annotation(&record, &ann))
        .for_each(|record| {
            fasta_writer
                .write_record(&record)
                .expect("Could not write record.")
        });
}
