use std;
use std::collections::HashMap;

#[derive(Debug)]
pub struct GenomesAndContigs {
    pub genomes: Vec<String>,
    pub contig_to_genome: Vec<String>
}


impl GenomesAndContigs {
    pub fn new() -> GenomesAndContigs {
        GenomesAndContigs {
            genomes: vec!(),
            contig_to_genome: vec!()
        }
    }

    pub fn establish_genome(&mut self, genome_name: String) -> usize {
        let index = self.genomes.len();
        self.genomes.push(genome_name);
        return index
    }

    pub fn insert(&mut self, contig_name: String) {
        self.contig_to_genome.push(contig_name);
    }
}
/// Finds the first occurence of element in a slice
fn find_first<T>(slice: &[T], element: T) -> Result<usize, &'static str>
    where T: std::cmp::PartialEq<T> {

    let mut index: usize = 0;
    for el in slice {
        if *el == element {
            return Ok(index)
            //let res: Result<usize, None> = Ok(index)
        }
        index += 1;
    }
    return Err("Element not found in slice")
}
