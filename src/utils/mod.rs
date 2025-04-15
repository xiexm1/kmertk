pub mod fasta;
pub mod kmer;

// 重导出常用功能
pub use fasta::{parse_fasta, fasta_to_map};
pub use kmer::{read_kmers, calculate_coverage, levenshtein}; 