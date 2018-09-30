// pub mod contig;
// pub mod genome;
// pub mod mosdepth_genome_coverage_estimators;
pub mod genomes_and_contigs;
// pub mod bam_generator;
// pub mod filter;
// pub mod external_command_checker;
// pub mod bwa_index_maintenance;

extern crate bio;
extern crate log;

extern crate rust_htslib;
extern crate env_logger;
extern crate nix;

// use bio::io::fasta::*;
use std::io;
use std::str;
use std::path::Path;
use genomes_and_contigs::GenomesAndContigs;

use std::cmp::min;
use std::collections;
use std::convert::AsRef;
use std::fs;
use std::io::prelude::*;


/// A FASTA record.
#[derive(Default, Clone, Debug)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
}

impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
        }
    }

    /// Create a new Fasta record from given attributes.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: String) -> Self {
        let desc = match desc {
            Some(desc) => Some(desc.to_owned()),
            _ => None,
        };
        Record {
            id: id.to_owned(),
            desc,
            seq: seq,
        }
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty()
    }

    /// Check validity of Fasta record.
    pub fn check(&self) -> Result<(), &str> {
        if self.id().is_empty() {
            return Err("Expecting id for Fasta record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }

        Ok(())
    }

    /// Return the id of the record.
    pub fn id(&self) -> &str {
        self.id.as_ref()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        match self.desc.as_ref() {
            Some(desc) => Some(&desc),
            None => None,
        }
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> String {
        self.seq.clone()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
    }
}

/// An iterator over the records of a Fasta file.
pub struct Records<R: io::Read> {
    reader: Reader<R>,
    error_has_occured: bool,
}

impl<R: io::Read> Iterator for Records<R> {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        if self.error_has_occured {
            None
        } else {
            let mut record = Record::new();
            match self.reader.read(&mut record) {
                Ok(()) if record.is_empty() => None,
                Ok(()) => Some(Ok(record)),
                Err(err) => {
                    self.error_has_occured = true;
                    Some(Err(err))
                }
            }
        }
    }
}

/// A FASTA reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    reader: io::BufReader<R>,
    line: String,
}

impl Reader<fs::File> {
    /// Read FASTA from given file path.
    pub fn from_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::open(path).map(Reader::new)
    }
}

impl<R: io::Read> Reader<R> {
    /// Create a new Fasta reader given an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let reader = Reader::new(fasta_file);
    /// # }
    /// ```
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }

    /// Read next FASTA record into the given `Record`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # use bio::io::fasta::Record;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// # let mut reader = Reader::new(fasta_file);
    /// let mut record = Record::new();
    /// reader.read(&mut record);
    ///
    /// assert_eq!(record.id(), "id");
    /// assert_eq!(record.desc().unwrap(), "desc");
    /// assert_eq!(record.seq().to_vec(), b"AAAA");
    /// # }
    /// ```
    pub fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if self.line.is_empty() {
            try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() {
                return Ok(());
            }
        }

        if !self.line.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }
        record.id = self.line[1..]
            .trim_right()
            .splitn(2, ' ')
            .nth(0)
            .map(|s| s.to_owned())
            .unwrap();
        record.desc = self.line[1..]
            .trim_right()
            .splitn(2, ' ')
            .nth(1)
            .map(|s| s.to_owned());
        loop {
            self.line.clear();
            try!(self.reader.read_line(&mut self.line));
            if self.line.is_empty() || self.line.starts_with('>') {
                break;
            }
            record.seq.push_str(&self.line);
        }

        Ok(())
    }

    /// Return an iterator over the records of this Fasta file.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # use bio::io::fasta::Record;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// # let reader = Reader::new(fasta_file);
    /// for record in reader.records() {
    ///     let record = record.unwrap();
    ///     assert_eq!(record.id(), "id");
    ///     assert_eq!(record.desc().unwrap(), "desc");
    ///     assert_eq!(record.seq().to_vec(), b"AAAA");
    /// }
    /// # }
    /// ```
    pub fn records(self) -> Records<R> {
        Records {
            reader: self,
            error_has_occured: false,
        }
    }
}



pub fn read_genome_fasta_files(fasta_file_paths: &[&str]) -> GenomesAndContigs {
    let mut contig_to_genome = GenomesAndContigs::new();
    let mut reader;
    for file in fasta_file_paths {
        let path = Path::new(*file);
        reader = Reader::from_file(path);
        let records = reader.unwrap().records();
        for record in records {
            let genome = record.unwrap();
            contig_to_genome.genomes.push(genome.id().to_string());
            contig_to_genome.contig_to_genome.push(genome.seq().to_string());
            
        }
    }
    return contig_to_genome;
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contig_to_genome(){
        let mut contig_to_genome = GenomesAndContigs::new();
        let genome = String::from("genome0");
        let index = contig_to_genome.establish_genome(genome);
        contig_to_genome.insert(String::from("contig1"), index);
        assert_eq!(
            String::from("genome0"),
            *(contig_to_genome.genome_of_contig(&String::from("contig1")).unwrap()));
    }

    #[test]
    fn test_read_genome_fasta_files_one_genome(){
        let contig_to_genome = read_genome_fasta_files(&vec!["tests/data/genome1.fna"]);
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq1")).unwrap());
        assert_eq!(String::from("genome1"), *contig_to_genome.genome_of_contig(&String::from("seq2")).unwrap());
    }
}
