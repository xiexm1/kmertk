use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use clap::Args;
use rayon::prelude::*;
use std::time::Instant;

use crate::utils::{read_kmers, parse_fasta, fasta_to_map};

#[derive(Args)]
pub struct CollectionArgs {
    /// 输入方式: 'seq' (序列文件构建矩阵), 'sam' (SAM文件), 'bowtie' (Bowtie文件)
    #[arg(short, long, default_value = "seq")]
    pub input_mode: String,

    /// 序列文件路径 (FASTA格式)
    #[arg(short = 'q', long)]
    pub sequence_file: Option<String>,

    /// kmer文件路径
    #[arg(short, long)]
    pub kmer_file: Option<String>,

    /// SAM或Bowtie格式的比对文件路径
    #[arg(short = 'a', long)]
    pub alignment_file: Option<String>,

    /// 输出文件路径
    #[arg(short, long)]
    pub output_file: String,

    /// 是否输出稀疏矩阵格式 (默认: false)
    #[arg(short = 'p', long)]
    pub sparse: bool,

    /// 最大允许的错配数量 (默认: 0)
    #[arg(short, long, default_value = "0")]
    pub max_mismatches: usize,

    /// 是否过滤CRISPR kmers (默认: false)
    #[arg(short = 'c', long, default_value = "false")]
    pub filter_crispr: bool,
}

// 稀疏矩阵结构
struct SparseMatrix {
    genes: Vec<String>,
    kmers_count: usize,
    rows: Vec<HashMap<usize, u8>>,
}

impl SparseMatrix {
    fn new(genes: &[String], kmers_count: usize) -> Self {
        let rows = vec![HashMap::new(); genes.len()];
        SparseMatrix { 
            genes: genes.to_vec(), 
            kmers_count, 
            rows 
        }
    }

    fn set(&mut self, row: usize, col: usize, value: u8) {
        if row < self.rows.len() {
            self.rows[row].insert(col, value);
        }
    }

    fn get(&self, row: usize, col: usize) -> u8 {
        if row < self.rows.len() {
            *self.rows[row].get(&col).unwrap_or(&0)
        } else {
            0
        }
    }
}

// SAM比对结构
struct SamAlignment {
    query_id: String,
    flag: u16,
    ref_id: String,
    position: usize,
    cigar: String,
    sequence: String,
    nm_value: usize,  // NM:i:值，表示错配数
}

impl SamAlignment {
    // 计算插入删除数量
    fn count_indels(&self) -> usize {
        let mut count = 0;
        let mut num_str = String::new();
        
        for c in self.cigar.chars() {
            if c.is_ascii_digit() {
                num_str.push(c);
            } else if c == 'I' || c == 'D' {
                if let Ok(num) = num_str.parse::<usize>() {
                    count += num;
                }
                num_str.clear();
            } else {
                num_str.clear();
            }
        }
        
        count
    }
    
    // 检查是否通过过滤条件
    fn passes_filter(&self, max_mismatches: usize) -> bool {
        // NM值包含所有错配和indels
        self.nm_value <= max_mismatches
    }
}

// Bowtie比对结构
struct BowtieAlignment {
    query_id: String,
    strand: char,
    ref_id: String,
    position: usize,
    sequence: String,
    mismatches: Vec<(usize, char, char)>, // 位置，参考碱基，查询碱基
}

// 近似匹配函数，允许最多max_mismatches个错配
fn find_approximate_match(text: &[u8], pattern: &[u8], max_mismatches: usize) -> bool {
    if pattern.len() > text.len() {
        return false;
    }
    
    // 遍历所有可能的起始位置
    for start in 0..=(text.len() - pattern.len()) {
        let mut mismatches = 0;
        
        // 检查当前起始位置的模式匹配
        for i in 0..pattern.len() {
            if text[start + i] != pattern[i] {
                mismatches += 1;
                if mismatches > max_mismatches {
                    break;
                }
            }
        }
        
        if mismatches <= max_mismatches {
            return true;
        }
    }
    
    false
}

// CRISPR匹配函数，检查k-mer是否为CRISPR相关序列
fn find_crispr_match(text: &[u8], pattern: &[u8]) -> bool {
    // CRISPR相关模式通常有重复模式和特定的间隔区
    // 这里简化实现，仅检查是否包含常见的CRISPR PAM序列
    let pam_sequences = [b"NGG", b"CCN", b"TTT"];
    
    for pam in &pam_sequences {
        for i in 0..=(pattern.len().saturating_sub(pam.len())) {
            let mut matches = true;
            for (j, &p) in pam.iter().enumerate() {
                if p != b'N' && pattern[i + j] != p {
                    matches = false;
                    break;
                }
            }
            if matches {
                return true;
            }
        }
    }
    
    false
}

// 过滤CRISPR相关k-mer
fn filter_crispr_kmers(kmers: &[String]) -> Vec<String> {
    kmers.par_iter()
        .filter(|kmer| !find_crispr_match(kmer.as_bytes(), kmer.as_bytes()))
        .cloned()
        .collect()
}

// 基于序列和k-mer构建矩阵
fn build_kmer_matrix(
    sequences: &HashMap<String, String>,
    kmers: &[String],
    max_mismatches: usize,
) -> SparseMatrix {
    let gene_ids: Vec<String> = sequences.keys().cloned().collect();
    let mut matrix = SparseMatrix::new(&gene_ids, kmers.len());
    
    println!("构建k-mer匹配矩阵: {} 个序列 x {} 个k-mer", sequences.len(), kmers.len());
    let start = Instant::now();
    
    // 并行处理每个序列
    let results: Vec<(usize, Vec<(usize, u8)>)> = gene_ids.par_iter().enumerate()
        .map(|(gene_idx, gene_id)| {
            let seq = sequences.get(gene_id).unwrap();
            let seq_bytes = seq.as_bytes();
            
            // 检查每个k-mer是否匹配当前序列
            let matches: Vec<(usize, u8)> = kmers.par_iter().enumerate()
                .filter_map(|(kmer_idx, kmer)| {
                    let kmer_bytes = kmer.as_bytes();
                    if find_approximate_match(seq_bytes, kmer_bytes, max_mismatches) {
                        Some((kmer_idx, 1))
                    } else {
                        None
                    }
                })
                .collect();
            
            (gene_idx, matches)
        })
        .collect();
    
    // 将结果填入矩阵
    for (gene_idx, matches) in results {
        for (kmer_idx, value) in matches {
            matrix.set(gene_idx, kmer_idx, value);
        }
    }
    
    println!("矩阵构建完成，耗时: {:?}", start.elapsed());
    matrix
}

// 从SAM文件构建矩阵
fn read_sam<P: AsRef<Path>>(sam_file: P) -> io::Result<Vec<SamAlignment>> {
    let file = File::open(sam_file)?;
    let reader = BufReader::new(file);
    let mut alignments = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        
        // 跳过头部注释行
        if line.starts_with('@') {
            continue;
        }
        
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue; // 跳过格式不正确的行
        }
        
        // 解析必要字段
        let query_id = fields[0].to_string();
        let flag = fields[1].parse::<u16>().unwrap_or(0);
        let ref_id = fields[2].to_string();
        let position = fields[3].parse::<usize>().unwrap_or(0);
        let cigar = fields[5].to_string();
        let sequence = fields[9].to_string();
        
        // 查找NM标签（错配数）
        let mut nm_value = 0;
        for i in 11..fields.len() {
            if fields[i].starts_with("NM:i:") {
                if let Some(val) = fields[i].strip_prefix("NM:i:") {
                    nm_value = val.parse::<usize>().unwrap_or(0);
                    break;
                }
            }
        }
        
        alignments.push(SamAlignment {
            query_id,
            flag,
            ref_id,
            position,
            cigar,
            sequence,
            nm_value,
        });
    }
    
    Ok(alignments)
}

// 从Bowtie文件构建矩阵
fn read_bowtie<P: AsRef<Path>>(bowtie_file: P) -> io::Result<Vec<BowtieAlignment>> {
    let file = File::open(bowtie_file)?;
    let reader = BufReader::new(file);
    let mut alignments = Vec::new();
    
    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        
        if fields.len() < 7 {
            continue; // 跳过格式不正确的行
        }
        
        let query_id = fields[0].to_string();
        let strand = if fields[1] == "+" { '+' } else { '-' };
        let ref_id = fields[2].to_string();
        let position = fields[3].parse::<usize>().unwrap_or(0);
        let sequence = fields[4].to_string();
        
        // 解析错配信息
        let mismatch_str = fields[7..].join(",");
        let mut mismatches = Vec::new();
        
        for m in mismatch_str.split(',') {
            let m_parts: Vec<&str> = m.split(':').collect();
            if m_parts.len() == 3 {
                if let Ok(pos) = m_parts[0].parse::<usize>() {
                    let ref_base = m_parts[1].chars().next().unwrap_or('N');
                    let query_base = m_parts[2].chars().next().unwrap_or('N');
                    mismatches.push((pos, ref_base, query_base));
                }
            }
        }
        
        alignments.push(BowtieAlignment {
            query_id,
            strand,
            ref_id,
            position,
            sequence,
            mismatches,
        });
    }
    
    Ok(alignments)
}

// 从SAM比对构建矩阵
fn build_sam_matrix(alignments: &[SamAlignment], max_mismatches: usize) -> (SparseMatrix, Vec<String>) {
    // 收集所有唯一的参考ID和查询ID
    let mut ref_ids = HashSet::new();
    let mut query_ids = Vec::new();
    
    for aln in alignments {
        if aln.passes_filter(max_mismatches) {
            ref_ids.insert(aln.ref_id.clone());
            if !query_ids.contains(&aln.query_id) {
                query_ids.push(aln.query_id.clone());
            }
        }
    }
    
    let ref_id_list: Vec<String> = ref_ids.into_iter().collect();
    let mut ref_id_map: HashMap<String, usize> = HashMap::new();
    for (i, id) in ref_id_list.iter().enumerate() {
        ref_id_map.insert(id.clone(), i);
    }
    
    // 创建矩阵
    let mut matrix = SparseMatrix::new(&ref_id_list, query_ids.len());
    
    // 填充矩阵
    for aln in alignments {
        if aln.passes_filter(max_mismatches) {
            let ref_idx = match ref_id_map.get(&aln.ref_id) {
                Some(&idx) => idx,
                None => continue,
            };
            
            let query_idx = match query_ids.iter().position(|id| id == &aln.query_id) {
                Some(idx) => idx,
                None => continue,
            };
            
            matrix.set(ref_idx, query_idx, 1);
        }
    }
    
    (matrix, query_ids)
}

// 从Bowtie比对构建矩阵
fn build_bowtie_matrix(alignments: &[BowtieAlignment]) -> (SparseMatrix, Vec<String>) {
    // 收集所有唯一的参考ID和查询ID
    let mut ref_ids = HashSet::new();
    let mut query_ids = Vec::new();
    
    for aln in alignments {
        ref_ids.insert(aln.ref_id.clone());
        if !query_ids.contains(&aln.query_id) {
            query_ids.push(aln.query_id.clone());
        }
    }
    
    let ref_id_list: Vec<String> = ref_ids.into_iter().collect();
    let mut ref_id_map: HashMap<String, usize> = HashMap::new();
    for (i, id) in ref_id_list.iter().enumerate() {
        ref_id_map.insert(id.clone(), i);
    }
    
    // 创建矩阵
    let mut matrix = SparseMatrix::new(&ref_id_list, query_ids.len());
    
    // 填充矩阵
    for aln in alignments {
        let ref_idx = match ref_id_map.get(&aln.ref_id) {
            Some(&idx) => idx,
            None => continue,
        };
        
        let query_idx = match query_ids.iter().position(|id| id == &aln.query_id) {
            Some(idx) => idx,
            None => continue,
        };
        
        matrix.set(ref_idx, query_idx, 1);
    }
    
    (matrix, query_ids)
}

// 输出完整矩阵
fn output_full_matrix(
    matrix: &SparseMatrix, 
    col_ids: &[String], 
    output_file: &str, 
    use_original_ids: bool
) -> io::Result<()> {
    let file = File::create(output_file)?;
    let mut writer = io::BufWriter::new(file);
    
    // 写入表头（列ID）
    if use_original_ids {
        write!(writer, "Gene_ID")?;
        for col_id in col_ids {
            write!(writer, ",{}", col_id)?;
        }
        writeln!(writer)?;
    } else {
        write!(writer, "Gene_ID")?;
        for col in 0..matrix.kmers_count {
            write!(writer, ",Col_{}", col)?;
        }
        writeln!(writer)?;
    }
    
    // 写入矩阵数据
    for (row, gene_id) in matrix.genes.iter().enumerate() {
        write!(writer, "{}", gene_id)?;
        
        for col in 0..matrix.kmers_count {
            let value = matrix.get(row, col);
            write!(writer, ",{}", value)?;
        }
        
        writeln!(writer)?;
    }
    
    Ok(())
}

// 输出稀疏矩阵
fn output_sparse_matrix(
    matrix: &SparseMatrix, 
    col_ids: &[String], 
    output_file: &str, 
    use_original_ids: bool
) -> io::Result<()> {
    let file = File::create(output_file)?;
    let mut writer = io::BufWriter::new(file);
    
    // 写入头部信息
    writeln!(writer, "# 稀疏矩阵格式: 每行包含 gene_id,kmer_id,value")?;
    writeln!(writer, "# 稀疏矩阵: {} 行 x {} 列", matrix.genes.len(), matrix.kmers_count)?;
    
    // 写入非零元素
    for (row, gene_id) in matrix.genes.iter().enumerate() {
        for col in 0..matrix.kmers_count {
            let value = matrix.get(row, col);
            if value != 0 {
                let col_id = if use_original_ids && col < col_ids.len() {
                    &col_ids[col]
                } else {
                    &format!("Col_{}", col)
                };
                
                writeln!(writer, "{},{},{}", gene_id, col_id, value)?;
            }
        }
    }
    
    Ok(())
}

pub fn run_collection(args: &CollectionArgs) -> io::Result<()> {
    match args.input_mode.as_str() {
        "seq" => {
            // 检查必要的参数
            let sequence_file = match &args.sequence_file {
                Some(file) => file,
                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "序列模式下需要提供序列文件路径 (--sequence-file)")),
            };
            
            let kmer_file = match &args.kmer_file {
                Some(file) => file,
                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "序列模式下需要提供k-mer文件路径 (--kmer-file)")),
            };
            
            // 读取序列和k-mer
            println!("读取序列文件: {}", sequence_file);
            let sequences = fasta_to_map(sequence_file)?;
            
            println!("读取k-mer文件: {}", kmer_file);
            let mut kmers = read_kmers(kmer_file)?;
            
            // 如果需要，过滤CRISPR相关k-mer
            if args.filter_crispr {
                println!("过滤前k-mer数量: {}", kmers.len());
                kmers = filter_crispr_kmers(&kmers);
                println!("过滤后k-mer数量: {}", kmers.len());
            }
            
            // 构建矩阵
            let matrix = build_kmer_matrix(&sequences, &kmers, args.max_mismatches);
            
            // 输出矩阵
            if args.sparse {
                println!("输出稀疏矩阵到: {}", args.output_file);
                output_sparse_matrix(&matrix, &kmers, &args.output_file, true)?;
            } else {
                println!("输出完整矩阵到: {}", args.output_file);
                output_full_matrix(&matrix, &kmers, &args.output_file, true)?;
            }
        },
        "sam" => {
            // 检查必要的参数
            let alignment_file = match &args.alignment_file {
                Some(file) => file,
                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "SAM模式下需要提供比对文件路径 (--alignment-file)")),
            };
            
            // 读取SAM文件
            println!("读取SAM文件: {}", alignment_file);
            let alignments = read_sam(alignment_file)?;
            println!("读取到 {} 条比对记录", alignments.len());
            
            // 构建矩阵
            let (matrix, query_ids) = build_sam_matrix(&alignments, args.max_mismatches);
            
            // 输出矩阵
            if args.sparse {
                println!("输出稀疏矩阵到: {}", args.output_file);
                output_sparse_matrix(&matrix, &query_ids, &args.output_file, true)?;
            } else {
                println!("输出完整矩阵到: {}", args.output_file);
                output_full_matrix(&matrix, &query_ids, &args.output_file, true)?;
            }
        },
        "bowtie" => {
            // 检查必要的参数
            let alignment_file = match &args.alignment_file {
                Some(file) => file,
                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Bowtie模式下需要提供比对文件路径 (--alignment-file)")),
            };
            
            // 读取Bowtie文件
            println!("读取Bowtie文件: {}", alignment_file);
            let alignments = read_bowtie(alignment_file)?;
            println!("读取到 {} 条比对记录", alignments.len());
            
            // 构建矩阵
            let (matrix, query_ids) = build_bowtie_matrix(&alignments);
            
            // 输出矩阵
            if args.sparse {
                println!("输出稀疏矩阵到: {}", args.output_file);
                output_sparse_matrix(&matrix, &query_ids, &args.output_file, true)?;
            } else {
                println!("输出完整矩阵到: {}", args.output_file);
                output_full_matrix(&matrix, &query_ids, &args.output_file, true)?;
            }
        },
        _ => {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                       "无效的输入模式，支持的模式: 'seq', 'sam', 'bowtie'"));
        }
    }
    
    Ok(())
} 