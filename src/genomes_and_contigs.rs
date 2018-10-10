
use std;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::str;

#[derive(Debug)]
pub struct GenomesAndContigs {
    pub genomes: Vec<String>,
    pub contig_to_genome: Vec<String>,
    pub kmers: Vec<HashMap<String, i32>>
}



impl GenomesAndContigs {
    pub fn new() -> GenomesAndContigs {
        GenomesAndContigs {
            genomes: vec!(),
            contig_to_genome: vec!(),
            kmers: vec!()
        }
    }

    pub fn establish_genome(&mut self, genome_name: String) -> usize {
        let index = self.genomes.len();
        self.genomes.push(genome_name);
        return index
    }

    pub fn establish_kmers(&mut self, kmer_size: usize){

        for (i, genome) in self.contig_to_genome.clone().iter().enumerate(){
            let mut start = 0;
            let mut stop = start + kmer_size;
            let mut full_gen = Vec::new();
            let mut gen = genome.lines();
            for line in gen{
                full_gen.push(line);
            }
            let mut kmer_hashmap = HashMap::new();
            let mut joined_genome = full_gen.join("");
            while stop <= joined_genome.len() as usize{
                let kmer = &joined_genome.as_bytes()[start..stop];
                let count = kmer_hashmap.entry(str::from_utf8(kmer).unwrap().to_string()).or_insert(0);
                *count += 1;
//                if self.kmers.contains_key(str::from_utf8(kmer).unwrap()){
//                    *self.kmers[i].get_mut(str::from_utf8(kmer).unwrap()).unwrap() += 1;
//                }else{
//                    self.kmers[i].insert(str::from_utf8(kmer).unwrap().to_string(), 1);
//                }
                start += 1;
                stop += 1;
            }
            self.kmers.push(kmer_hashmap);

        }
    }

    pub fn insert(&mut self, contig_name: String) {
        self.contig_to_genome.push(contig_name);
    }

    pub fn get_kmers(&mut self){
        print!("K-Mer");
        for genome in self.genomes.clone(){
            print!("\t{}", genome.split("/").collect::<Vec<_>>().last().unwrap());
        }
        print!("\n");
        let mut full_kmer_map: HashMap<String, Vec<i32>> = HashMap::new();
        for (index, kmer_map) in self.kmers.iter().enumerate(){
            for (key, value) in kmer_map.iter() {
                let value_vec = full_kmer_map.entry(key.clone()).or_insert(vec![0; self.genomes.len()]);
                value_vec[index] = *value;
            }
        }
        for (key, value) in full_kmer_map{
            print!("{}", key);
            for v in value{
                print!("\t{}", v);
            }
            print!("\n")
        }
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
