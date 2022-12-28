use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Format, Header, Read, Reader, Writer};
use std::env;
use std::process;
use std::str;

fn filter_cigars(input_bam_filename: &str, output_bam_filename: &str) {
    let mut input_bam = Reader::from_path(input_bam_filename).unwrap();
    let header = Header::from_template(input_bam.header());
    let mut output_bam = Writer::from_path(output_bam_filename, &header, Format::BAM).unwrap();

    input_bam
        .set_threads(2)
        .expect("Failed to set number of BAM reading threads to 2.");
    output_bam
        .set_threads(5)
        .expect("Failed to set number of BAM writing threads to 4.");

    for r in input_bam.records() {
        let record = r.unwrap();
        let ok = check_cigar(&record);
        if ok {
            output_bam.write(&record).unwrap();
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 3 {
        eprintln!("Usage: filter_cigars input.bam output.bam");
        process::exit(1);
    }

    let input_bam_filename = &args[1];
    let output_bam_filename = &args[2];

    filter_cigars(input_bam_filename, output_bam_filename);
}

fn check_cigar(read: &bam::Record) -> bool {
    // These variables keep track of our current position
    // within the pileup track and the read sequence, as we
    // read through the CIGAR string.
    let mut pileup_idx = 0;
    let mut ok = true;

    for s in read.cigar().iter() {
        match *s {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                pileup_idx += len;
            }
            Cigar::Ins(_) => {
                if pileup_idx == 0 {
                    ok = false;
                }
            }
            Cigar::Del(_) => {
                if pileup_idx == 0 {
                    ok = false;
                }
            }
            _ => {
                ok = false;
            }
        }
    }
    ok
}
