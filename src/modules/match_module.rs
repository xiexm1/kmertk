use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::Path;
use clap::Args;
use aho_corasick::{AhoCorasick, PatternID};
use rayon::prelude::*;
use flate2::read::GzDecoder;
use bio::io::fastq;
use std::io::BufWriter;

use crate::utils::{read_kmers, parse_fasta};

#[derive(Args)]
pub struct MatchArgs {
    /// K-mer文件路径
    #[arg(short, long)]
    pub kmer_file: String,

    /// 输入文件路径 (支持 FASTA, FASTQ, 和压缩文件 .gz)
    #[arg(short, long)]
    pub input_file: String,

    /// 输出文件路径 (匹配结果)
    #[arg(short, long)]
    pub output_file: String,

    /// 是否提取匹配序列 (默认: false)
    #[arg(short, long, default_value = "false")]
    pub extract_sequences: bool,

    /// 提取序列输出文件路径 (当 extract_sequences 为 true 时使用)
    #[arg(short = 'x', long)]
    pub extract_file: Option<String>,
    
    /// 输出格式: "fasta" 或 "fastq" (默认: 与输入格式相同)
    #[arg(short = 'f', long, default_value = "same")]
    pub output_format: String,
}

// 支持的文件类型
enum FileType {
    Fasta,
    Fastq,
}

pub fn run_match(args: &MatchArgs) -> io::Result<()> {
    // 检查输入文件是否存在
    if !Path::new(&args.kmer_file).exists() {
        return Err(io::Error::new(io::ErrorKind::NotFound, format!("K-mer文件不存在: {}", args.kmer_file)));
    }
    if !Path::new(&args.input_file).exists() {
        return Err(io::Error::new(io::ErrorKind::NotFound, format!("输入文件不存在: {}", args.input_file)));
    }

    // 如果需要提取序列但没有指定提取文件路径，返回错误
    if args.extract_sequences && args.extract_file.is_none() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, 
            "启用提取序列功能时，必须指定提取文件路径 (--extract_file)"));
    }

    // 读取k-mers
    let kmers = read_kmers(&args.kmer_file)?;

    // 构建 Aho-Corasick 自动机
    let ac = match AhoCorasick::new(&kmers) {
        Ok(automaton) => automaton,
        Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, 
                                           format!("Aho-Corasick 构建失败: {}", e))),
    };

    // 确定输入文件类型
    let file_type = if args.input_file.ends_with(".fa") || args.input_file.ends_with(".fasta") 
        || args.input_file.ends_with(".fa.gz") || args.input_file.ends_with(".fasta.gz") {
        FileType::Fasta
    } else {
        FileType::Fastq
    };

    // 存储结果
    let mut matches: HashMap<String, HashSet<PatternID>> = HashMap::new();
    
    // 用于存储要提取的序列
    let mut extracted_sequences: Vec<(String, String, Option<String>)> = Vec::new(); // (id, seq, qual)

    // 处理输入文件
    match file_type {
        FileType::Fasta => process_fasta(&args.input_file, &ac, &kmers, &mut matches, 
                                         args.extract_sequences, &mut extracted_sequences)?,
        FileType::Fastq => process_fastq(&args.input_file, &ac, &kmers, &mut matches, 
                                         args.extract_sequences, &mut extracted_sequences)?,
    };

    // 写入匹配结果到文件
    let mut output = File::create(&args.output_file)?;
    for (seq_id, kmer_indices) in &matches {
        writeln!(output, "{}:", seq_id)?;
        for &idx in kmer_indices {
            writeln!(output, "  - {}", kmers[idx.as_usize()])?;
        }
        writeln!(output)?;
    }

    // 如果需要，提取匹配的序列到新文件
    if args.extract_sequences && !extracted_sequences.is_empty() {
        let extract_path = args.extract_file.as_ref().unwrap().clone();
        
        // 确定输出格式
        let output_format = match args.output_format.as_str() {
            "fasta" => FileType::Fasta,
            "fastq" => FileType::Fastq,
            "same" => file_type,
            _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                             "输出格式无效。有效选项: 'fasta', 'fastq', 'same'")),
        };
        
        // 写入提取的序列
        match output_format {
            FileType::Fasta => write_fasta(&extract_path, &extracted_sequences)?,
            FileType::Fastq => write_fastq(&extract_path, &extracted_sequences)?,
        }
    }

    println!("处理完成！结果已保存到: {}", args.output_file);
    if args.extract_sequences {
        println!("匹配的序列已提取到: {}", args.extract_file.as_ref().unwrap());
    }
    
    Ok(())
}

// 处理FASTA文件（包括gzip压缩文件）
fn process_fasta(
    input_file: &str, 
    ac: &AhoCorasick, 
    _kmers: &[String], 
    matches: &mut HashMap<String, HashSet<PatternID>>,
    extract: bool,
    extracted: &mut Vec<(String, String, Option<String>)>
) -> io::Result<()> {
    let sequences = parse_fasta(input_file)?;

    // 并行处理序列
    let parallel_matches: Vec<(String, HashSet<PatternID>, String)> = sequences.par_iter()
        .map(|(id, seq)| {
            let mut seq_matches = HashSet::new();
            // 使用find_iter查找匹配
            ac.find_iter(seq).for_each(|mat| {
                seq_matches.insert(mat.pattern());
            });
            (id.clone(), seq_matches, seq.clone())
        })
        .filter(|(_, seq_matches, _)| !seq_matches.is_empty())
        .collect();

    // 合并结果
    for (id, matches_set, seq) in parallel_matches {
        matches.insert(id.clone(), matches_set);
        if extract {
            extracted.push((id, seq, None));
        }
    }

    Ok(())
}

// 处理FASTQ文件（包括gzip压缩文件）
fn process_fastq(
    input_file: &str, 
    ac: &AhoCorasick, 
    _kmers: &[String],
    matches: &mut HashMap<String, HashSet<PatternID>>,
    extract: bool,
    extracted: &mut Vec<(String, String, Option<String>)>
) -> io::Result<()> {
    let reader: Box<dyn io::Read> = if input_file.ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(input_file)?))
    } else {
        Box::new(File::open(input_file)?)
    };
    
    let buf_reader = BufReader::new(reader);
    let fastq_reader = fastq::Reader::new(buf_reader);
    
    let mut sequences: Vec<(String, String, String)> = Vec::new();
    
    for record_result in fastq_reader.records() {
        // 处理fastq::Error到io::Error的转换
        let record = match record_result {
            Ok(rec) => rec,
            Err(e) => return Err(io::Error::new(io::ErrorKind::InvalidData, 
                                              format!("FASTQ解析错误: {}", e))),
        };
        
        let id = record.id().to_string();
        let seq = String::from_utf8_lossy(record.seq()).to_string();
        let qual = String::from_utf8_lossy(record.qual()).to_string();
        sequences.push((id, seq, qual));
    }
    
    // 并行处理序列
    let parallel_matches: Vec<(String, HashSet<PatternID>, String, String)> = sequences.par_iter()
        .map(|(id, seq, qual)| {
            let mut seq_matches = HashSet::new();
            // 使用find_iter查找匹配
            ac.find_iter(seq).for_each(|mat| {
                seq_matches.insert(mat.pattern());
            });
            (id.clone(), seq_matches, seq.clone(), qual.clone())
        })
        .filter(|(_, seq_matches, _, _)| !seq_matches.is_empty())
        .collect();

    // 合并结果
    for (id, matches_set, seq, qual) in parallel_matches {
        matches.insert(id.clone(), matches_set);
        if extract {
            extracted.push((id, seq, Some(qual)));
        }
    }

    Ok(())
}

// 写入FASTA格式文件
fn write_fasta(output_file: &str, sequences: &[(String, String, Option<String>)]) -> io::Result<()> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);
    
    for (id, seq, _) in sequences {
        writeln!(writer, ">{}", id)?;
        for chunk in seq.as_bytes().chunks(80) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap())?;
        }
    }
    
    Ok(())
}

// 写入FASTQ格式文件
fn write_fastq(output_file: &str, sequences: &[(String, String, Option<String>)]) -> io::Result<()> {
    let file = File::create(output_file)?;
    let mut writer = BufWriter::new(file);
    
    for (id, seq, qual_opt) in sequences {
        writeln!(writer, "@{}", id)?;
        writeln!(writer, "{}", seq)?;
        writeln!(writer, "+")?;
        
        if let Some(qual) = qual_opt {
            writeln!(writer, "{}", qual)?;
        } else {
            // 如果没有质量分数，生成默认的质量分数（全部为'I'，表示较高质量）
            let default_qual = "I".repeat(seq.len());
            writeln!(writer, "{}", default_qual)?;
        }
    }
    
    Ok(())
} 