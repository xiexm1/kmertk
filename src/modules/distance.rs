use std::io::{self, Write};
use std::fs::File;
use std::path::Path;
use clap::Args;
use rand::seq::SliceRandom;
use rayon::prelude::*;

use crate::utils::{parse_fasta, levenshtein};

#[derive(Args)]
pub struct DistanceArgs {
    /// 序列FASTA文件路径
    #[arg(short, long)]
    pub input_file: String,

    /// 要计算的序列对数量
    #[arg(short, long, default_value = "100")]
    pub pairs: usize,

    /// 输出文件路径
    #[arg(short, long)]
    pub output_file: String,
    
    /// 是否使用所有可能的序列对（而不是随机采样）
    #[arg(short, long, default_value = "false")]
    pub all_pairs: bool,
}

pub fn run_distance(args: &DistanceArgs) -> io::Result<()> {
    // 检查输入文件是否存在
    if !Path::new(&args.input_file).exists() {
        return Err(io::Error::new(io::ErrorKind::NotFound, 
                                 format!("输入文件不存在: {}", args.input_file)));
    }

    // 读取序列
    let sequences = parse_fasta(&args.input_file)?;
    
    if sequences.len() < 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidData, 
                                 "FASTA文件必须包含至少两个序列才能计算距离。"));
    }
    
    let results = if args.all_pairs {
        // 计算所有可能的序列对距离
        compute_all_pairs(&sequences)
    } else {
        // 随机选择n对序列计算距离
        compute_random_pairs(&sequences, args.pairs)
    };

    // 写入结果到输出文件
    let mut output = File::create(&args.output_file)?;
    for (id1, id2, distance) in results {
        writeln!(output, "{}\t{}\t{}", id1, id2, distance)?;
    }

    // 输出统计信息
    let pairs_str = if args.all_pairs { "所有".to_string() } else { args.pairs.to_string() };
    println!("已计算{}对序列的距离，结果保存至: {}", pairs_str, args.output_file);

    Ok(())
}

// 计算随机选取的n对序列的距离
fn compute_random_pairs(sequences: &[(String, String)], n: usize) -> Vec<(String, String, usize)> {
    let mut pairs = Vec::new();
    
    for _ in 0..n {
        let idx1 = rand::random::<usize>() % sequences.len();
        let mut idx2 = rand::random::<usize>() % sequences.len();
        
        // 确保不会比较相同的序列
        while idx1 == idx2 {
            idx2 = rand::random::<usize>() % sequences.len();
        }
        
        pairs.push((&sequences[idx1], &sequences[idx2]));
    }
    
    pairs.par_iter()
        .map(|(&(ref id1, ref seq1), &(ref id2, ref seq2))| {
            let distance = levenshtein(seq1, seq2);
            (id1.clone(), id2.clone(), distance)
        })
        .collect()
}

// 计算所有可能的序列对的距离（注意：对于大量序列会有O(n²)的复杂度）
fn compute_all_pairs(sequences: &[(String, String)]) -> Vec<(String, String, usize)> {
    let mut pairs = Vec::new();
    let n = sequences.len();
    
    // 生成所有可能的序列对组合（不重复）
    for i in 0..n {
        for j in (i+1)..n {
            pairs.push((&sequences[i], &sequences[j]));
        }
    }
    
    // 并行计算距离
    pairs.par_iter()
        .map(|(&(ref id1, ref seq1), &(ref id2, ref seq2))| {
            let distance = levenshtein(seq1, seq2);
            (id1.clone(), id2.clone(), distance)
        })
        .collect()
} 