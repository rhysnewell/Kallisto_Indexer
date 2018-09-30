extern crate kallisto_indexer;

extern crate bio;
extern crate clap;
use clap::*;
extern crate rust_htslib;
extern crate nix;
extern crate tempdir;
extern crate csv;

extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

use std::env;
use std::str;
use std::process;
use std::path::Path;
use std::error::Error;
use std::ffi::OsString;
use std::fs::File;

fn main() {

    let mut app = build_cli();
    let m = app.clone().get_matches();
    let genomes_and_contigs;
    if m.is_present("fasta-files"){
        let genome_fasta_files: Vec<&str> = m.values_of("fasta-files").unwrap().collect();
        genomes_and_contigs = kallisto_indexer::read_genome_fasta_files(&genome_fasta_files);
    } else if m.is_present("fasta-directory") {
        let file_path = m.value_of("fasta-directory").unwrap();
        let mut rdr = csv::Reader::from_path(file_path);
        let mut genome_fasta_files: Vec<String> = vec!();
        for line in rdr.unwrap().records() {
            let record = line.unwrap();
            let file = &record[1];
            let mut s = String::from(file);
            genome_fasta_files.push(s);
        }
        let mut strs: Vec<&str> = vec!();
        for f in &genome_fasta_files {
            strs.push(f);
        }
        genomes_and_contigs = kallisto_indexer::read_genome_fasta_files(&strs);
        let mut genomes_string = String::new();
        for (key, value) in genomes_and_contigs.contig_to_genome{
            genomes_string = [genomes_string, key].join("N");
        }

    }
}

fn build_cli() -> App<'static, 'static> {

    return App::new("kallisto_indexer")
        .author("Rhys J. P. Newell <r.newell near uq.edu.au>")
        .about("Mapping coverage analysis for metagenomics")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .help("")
        .arg(Arg::with_name("fasta-files")
                .short("f")
                .long("fasta-files")
                .conflicts_with("fasta-directory")
                .multiple(true)
                .takes_value(true))
        .arg(Arg::with_name("threads")
                .short("-t")
                .long("threads")
                .default_value("1")
                .takes_value(true))
        .arg(Arg::with_name("fasta-directory")
                .short("d")
                .long("fasta-directory")
                .conflicts_with("fasta-files")
                .takes_value(true))
        .arg(Arg::with_name("verbose")
                .short("v")
                .long("verbose"))
        .arg(Arg::with_name("quiet")
                .short("q")
                .long("quiet"))
}