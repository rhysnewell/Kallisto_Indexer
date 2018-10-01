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
use std::io::{self, Write, Read};
use std::env;
use std::str;
use std::process::{self, Command};
use std::path::Path;
use std::error::Error;
use std::ffi::OsString;
use std::fs::{self, File};
use csv::ReaderBuilder;

fn main() {

    let mut app = build_cli();
    let m = app.clone().get_matches();
    let genomes_and_contigs;
    let mut genomes_string = String::new();
    if m.is_present("fasta-files"){
        let genome_fasta_files: Vec<&str> = m.values_of("fasta-files").unwrap().collect();
        genomes_and_contigs = kallisto_indexer::read_genome_fasta_files(&genome_fasta_files);
        for (i, value) in genomes_and_contigs.contig_to_genome.iter().enumerate(){
            genomes_string = [genomes_string, value.clone()].join(&format!(">{}\n", genomes_and_contigs.genomes[i]));
            genomes_string.push_str("\n")
        }
        run_kallisto(genomes_string, m.clone());
    } else if m.is_present("fasta-directory") {
        // println!("{}",format!("{:?}", m.value_of("fasta-directory").unwrap()));
        let mut file_path = File::open(m.value_of("fasta-directory").unwrap()).unwrap();
        // println!("{}",format!("{:?}", file_path));
        let mut contents = String::new();
        file_path.read_to_string(&mut contents);
        // println!("{}",format!("{:?}", contents));
        let mut rdr = csv::ReaderBuilder::new()
                                        .delimiter(b'\t')
                                        .has_headers(false)
                                        .from_reader(contents.as_bytes());
        let mut genome_fasta_files: Vec<String> = Vec::new();
        for line in rdr.records() {
            let record = line.unwrap();
            // println!("{}",format!("{:?}", record));
            let file = &record[1];
            // let mut s = String::from(file);
            genome_fasta_files.push(file.to_string());
        }
        let mut strs: Vec<&str> = vec!();
        for f in &genome_fasta_files {
            strs.push(f);
        }
        // info!("{:?}", &strs);
        genomes_and_contigs = kallisto_indexer::read_genome_fasta_files(&strs);
        for (i, value) in genomes_and_contigs.contig_to_genome.iter().enumerate(){
            genomes_string = [genomes_string, value.clone()].join(&format!(">{}\n", genomes_and_contigs.genomes[i]));
            genomes_string.push_str("\n")
        }
        run_kallisto(genomes_string, m.clone());

    }
    
}

fn run_kallisto(genomes_string: String, matches: ArgMatches) -> Result<()>{
        let m = matches;
        let mkdir = Command::new("mkdir")
                .arg(format!("{}", m.value_of("output").unwrap()))
                .output()
                .expect("Failed");
        println!("stderr mkdir: {}", String::from_utf8(mkdir.stderr).unwrap());
        let dir = Builder::new().tempdir_in(format!("{}", m.value_of("output").unwrap()))?;
        // info!("Temp Path {:?}", dir.unwrap().path());
        let temp_file_path = dir.path().join("temp_genome.fasta");
        let mut temp_file = File::create(temp_file_path.clone())?;
        let mut output;
        // println!("{}", &genomes_string);
        temp_file.write_all(genomes_string.as_bytes())?;
        let touch = Command::new("touch")
                .arg(format!("{}/genomes.idx", m.value_of("output").unwrap()))
                .output()
                .expect("Failed");
        println!("stderr touch: {}", String::from_utf8(touch.stderr).unwrap());
        if m.is_present("k-mer-size"){
            output = Command::new("kallisto")
                    .arg("index")
                    .arg(format!("--index={}/genomes.idx", m.value_of("output").unwrap()))
                    .arg(format!("--kmer-size={}", m.value_of("k-mer-size").unwrap()))
                    .arg(temp_file_path)
                    .output()
                    .expect("failed to execute process");
        }else{
            
            output = Command::new("kallisto")
                    .arg("index")
                    .arg(format!("--index={}/genomes.idx", m.value_of("output").unwrap()))
                    .arg(temp_file_path)
                    .output()
                    .expect("failed to execute process");
        }
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stdout: {}", String::from_utf8_lossy(&output.stderr));
        // drop(temp_file);
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
        .arg(Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .required(true))
        .arg(Arg::with_name("verbose")
                .short("v")
                .long("verbose"))
        .arg(Arg::with_name("quiet")
                .short("q")
                .long("quiet"))
}