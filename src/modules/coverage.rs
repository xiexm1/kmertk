use std::io::{self, Write};
use std::fs::File;
use std::path::Path;
use std::collections::HashSet;
use clap::Args;
use rayon::prelude::*;

use crate::utils::{parse_fasta, read_kmers, calculate_coverage};

#[derive(Args)]
pub struct CoverageArgs {
    /// 序列FASTA文件路径
    #[arg(short, long)]
    pub input_file: String,

    /// k-mer文件路径
    #[arg(short, long)]
    pub kmer_file: String,

    /// 输出文件路径
    #[arg(short, long)]
    pub output_file: String,
    
    /// 只显示覆盖率超过阈值的序列
    #[arg(short, long, default_value = "0.0")]
    pub threshold: f64,
}

pub fn run_coverage(args: &CoverageArgs) -> io::Result<()> {
    // 检查输入文件是否存在
    if !Path::new(&args.input_file).exists() {
        return Err(io::Error::new(io::ErrorKind::NotFound, 
                                 format!("输入文件不存在: {}", args.input_file)));
    }
    if !Path::new(&args.kmer_file).exists() {
        return Err(io::Error::new(io::ErrorKind::NotFound, 
                                 format!("k-mer文件不存在: {}", args.kmer_file)));
    }

    // 读取序列
    let sequences = parse_fasta(&args.input_file)?;
    
    // 读取k-mer并创建集合
    let kmers = read_kmers(&args.kmer_file)?;
    let kmer_set: HashSet<&str> = kmers.iter().map(|s| s.as_str()).collect();
    
    // 并行计算每个序列的覆盖率
    let results: Vec<(String, f64)> = sequences.par_iter()
        .map(|(id, seq)| {
            let coverage = calculate_coverage(seq, &kmer_set);
            (id.clone(), coverage)
        })
        .filter(|(_, coverage)| *coverage >= args.threshold)
        .collect();
    
    // 写入结果到输出文件
    let mut output = if args.output_file == "-" {
        // 如果输出文件为"-"，则输出到标准输出
        Box::new(io::stdout()) as Box<dyn Write>
    } else {
        Box::new(File::create(&args.output_file)?) as Box<dyn Write>
    };
    
    // 写入表头
    writeln!(output, "序列ID\t覆盖率\t覆盖率百分比")?;
    
    // 写入每个序列的覆盖率
    for (id, coverage) in &results {
        writeln!(output, "{}\t{:.4}\t{:.2}%", id, coverage, coverage * 100.0)?;
    }
    
    // 输出统计信息
    if args.output_file != "-" {
        println!("已计算{}个序列的覆盖率，结果保存至: {}", results.len(), args.output_file);
    }
    
    // 计算并显示汇总统计
    let avg_coverage: f64 = if results.is_empty() {
        0.0
    } else {
        results.iter().map(|(_, cov)| cov).sum::<f64>() / results.len() as f64
    };
    
    println!("覆盖率统计:");
    println!("  - 序列总数: {}", sequences.len());
    println!("  - 超过阈值的序列数: {}", results.len());
    println!("  - 平均覆盖率: {:.2}%", avg_coverage * 100.0);
    
    Ok(())
} 