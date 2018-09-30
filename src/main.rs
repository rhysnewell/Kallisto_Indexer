extern crate kallisto_indexer;

extern crate bio;
extern crate clap;
use clap::*;
extern crate rust_htslib;
extern crate nix;
extern crate tempfile;
extern crate csv;

#[macro_use]
extern crate log;
use log::LogLevelFilter;
extern crate env_logger;
use env_logger::LogBuilder;

use tempfile::{tempfile, tempdir, Builder};
use std::io::{self, Write};
use std::env;
use std::str;
use std::process::{self, Command};
use std::path::Path;
use std::error::Error;
use std::ffi::OsString;
use std::fs::{self, File};

fn main() {

    let mut app = build_cli();
    let m = app.clone().get_matches();
    let genomes_and_contigs;
    let mut genomes_string = String::new();
    if m.is_present("fasta-files"){
        let genome_fasta_files: Vec<&str> = m.values_of("fasta-files").unwrap().collect();
        genomes_and_contigs = kallisto_indexer::read_genome_fasta_files(&genome_fasta_files);
        for (key, value) in genomes_and_contigs.contig_to_genome{
            genomes_string = [genomes_string, key].join("N");
        }
        run_kallisto(genomes_string, m.clone());
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
        for (key, value) in genomes_and_contigs.contig_to_genome{
            genomes_string = [genomes_string, key].join("N");
        }
        run_kallisto(genomes_string, m.clone());

    }
    
}

fn run_kallisto(genomes_string: String, matches: ArgMatches) -> Result<()>{
        let m = matches;
        let dir = Builder::new().tempdir_in("./")?;
        // info!("Temp Path {:?}", dir.unwrap().path());
        let temp_file_path = dir.path().join("temp_genome.fasta");
        let mut temp_file = File::create(temp_file_path.clone())?;
        let mut output;
        println!("{}", &genomes_string);
        writeln!(temp_file, "{}", genomes_string)?;
        if m.is_present("k-mer-size"){
            output = Command::new("kallisto")
                    .arg("index")
                    .arg("-i")
                    .arg(temp_file_path)
                    .arg("-k")
                    .arg(m.value_of("k-mer-size").unwrap())
                    .output();
        }else{
            output = Command::new("kallisto")
                    .arg("index")
                    .arg("-i")
                    .arg(temp_file_path)
                    .output();
        }
        println!("stdout: {}", String::from_utf8_lossy(&output.unwrap().stdout));
        drop(temp_file);
        Ok(dir.close()?)
    }

fn build_cli() -> App<'static, 'static> {

    return App::new("kallisto_indexer")
        .author("Rhys J. P. Newell <r.newell near uq.edu.au>")
        .about("Mapping coverage analysis for metagenomics")
        .help("")
        .arg(Arg::with_name("fasta-files")
                .short("f")
                .long("fasta-files")
                .conflicts_with("fasta-directory")
                .multiple(true)
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("threads")
                .short("-t")
                .long("threads")
                .default_value("1")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("fasta-directory")
                .short("d")
                .long("fasta-directory")
                .conflicts_with("fasta-files")
                .takes_value(true))
        .arg(Arg::with_name("k-mer-size")
                .short("k")
                .long("k-mer-size")
                .takes_value(true)
                .required(false))
        .arg(Arg::with_name("verbose")
                .short("v")
                .long("verbose"))
        .arg(Arg::with_name("quiet")
                .short("q")
                .long("quiet"))
}