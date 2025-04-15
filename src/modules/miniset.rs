use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use clap::Args;
use rayon::prelude::*;
use bit_vec::BitVec;
use indicatif::{ProgressBar, ProgressStyle};
use std::time::Instant;

use crate::utils::read_kmers;

#[derive(Args)]
pub struct MiniSetArgs {
    /// 矩阵文件路径 (基因-kmer覆盖矩阵)
    #[arg(short, long)]
    pub matrix_file: String,

    /// 输出前缀 (用于生成多个输出文件)
    #[arg(short, long, default_value = "kmer_selection")]
    pub output_prefix: String,
    
    /// 目标覆盖率 (百分比)
    #[arg(short, long, default_value = "90.0")]
    pub target_coverage: f64,
    
    /// 是否使用稀疏矩阵格式
    #[arg(short, long, default_value = "false")]
    pub sparse: bool,
    
    /// 使用的线程数
    #[arg(short, long, default_value = "8")]
    pub threads: usize,
}

// 数据结构，使用BitVec提高内存效率和计算速度
struct Data {
    genes: Vec<String>,           // 基因名称
    kmers: Vec<String>,           // k-mer 名称
    matrix: Vec<BitVec>,          // 稀疏二进制矩阵（行：基因，列：k-mer）
    gene_hits: Vec<usize>,        // 记录每个基因被多少k-mer覆盖
}

pub fn run_miniset(args: &MiniSetArgs) -> io::Result<()> {
    // 检查输入文件是否存在
    if !Path::new(&args.matrix_file).exists() {
        return Err(io::Error::new(io::ErrorKind::NotFound, 
                                 format!("矩阵文件不存在: {}", args.matrix_file)));
    }
    
    // 设置线程数
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap_or_else(|e| eprintln!("设置线程池失败: {}", e));
    
    println!("使用 {} 个线程", rayon::current_num_threads());
    
    // 创建日志文件
    let log_file = format!("{}_selection.log", args.output_prefix);
    let mut logger = File::create(&log_file)?;
    
    let start_time = Instant::now();
    writeln!(logger, "开始时间: {:?}", start_time)?;
    
    // 根据输入格式加载数据
    let data = if args.sparse {
        load_sparse_matrix(&args.matrix_file, &mut logger, &args.output_prefix)?
    } else {
        load_matrix(&args.matrix_file, &mut logger, &args.output_prefix)?
    };
    
    println!("矩阵加载完成: {} 个基因, {} 个k-mer", data.genes.len(), data.kmers.len());
    
    // 执行贪心选择算法
    let target_counts = vec![50, 100, 200, 500, 1000];
    let results = greedy_kmer_selection(&data, &target_counts, &mut logger, &args.output_prefix)?;
    
    // 输出结果
    let stats_file = format!("{}_coverage.stats.txt", args.output_prefix);
    let mut stats_writer = File::create(&stats_file)?;
    writeln!(stats_writer, "kmer_count\tcovered_genes\tcoverage_percent\tdelta_coverage")?;
    
    println!("\n【覆盖结果统计】");
    println!("k-mer数量\t覆盖基因数\t覆盖率\t增长率");
    
    let mut last_coverage = 0.0;
    for (i, result) in results.iter().enumerate() {
        let (target, covered, _selected_kmers) = result;
        let coverage_percent = (*covered as f64 / data.genes.len() as f64) * 100.0;
        let delta_coverage = if i == 0 { coverage_percent } else { coverage_percent - last_coverage };
        
        // 写入统计文件
        writeln!(stats_writer, "{}\t{}\t{:.2}\t{:.2}", target, covered, coverage_percent, delta_coverage)?;
        
        println!("{}\t{}\t{:.2}%\t{:+.2}%", target, covered, coverage_percent, delta_coverage);
        
        last_coverage = coverage_percent;
    }
    
    // 保存最终选择的k-mer
    if let Some((_, _, selected_kmers)) = results.last() {
        let selected_file = format!("{}_selected_kmers.txt", args.output_prefix);
        let mut kmer_writer = File::create(&selected_file)?;
        
        for &kmer_idx in selected_kmers {
            writeln!(kmer_writer, "{}", data.kmers[kmer_idx])?;
        }
        
        println!("已保存选择的k-mer到: {}", selected_file);
    }
    
    let elapsed = start_time.elapsed();
    println!("总运行时间: {:?}", elapsed);
    writeln!(logger, "总运行时间: {:?}", elapsed)?;
    
    Ok(())
}

fn load_matrix(file_path: &str, logger: &mut File, output_prefix: &str) -> io::Result<Data> {
    let start = Instant::now();
    println!("开始读取矩阵文件...");
    writeln!(logger, "开始读取矩阵文件: {}", file_path)?;
    
    // 打开CSV文件
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    
    // 读取头行获取k-mer列表
    let header_line = match lines.next() {
        Some(Ok(line)) => line,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidData, "文件为空或无法读取头行")),
    };
    
    let headers: Vec<&str> = header_line.split(',').collect();
    let kmers = headers.iter().skip(1).map(|&s| s.to_string()).collect::<Vec<_>>();
    let num_kmers = kmers.len();
    
    println!("检测到 {} 个k-mer列", num_kmers);
    writeln!(logger, "检测到 {} 个k-mer列", num_kmers)?;
    
    let mut genes = Vec::new();
    let mut matrix = Vec::new();
    let mut gene_hits = Vec::new();
    let mut zero_count = 0;
    
    for (i, line_result) in lines.enumerate() {
        let line = match line_result {
            Ok(l) => l,
            Err(e) => {
                writeln!(logger, "警告: 第{}行读取错误: {}", i+2, e)?;
                continue;
            }
        };
        
        let values: Vec<&str> = line.split(',').collect();
        if values.len() <= 1 {
            writeln!(logger, "警告: 第{}行格式不正确: {}", i+2, line)?;
            continue;
        }
        
        let gene_id = values[0].to_string();
        genes.push(gene_id);
        
        let mut bitvec = BitVec::from_elem(num_kmers, false);
        let mut hit_count = 0;
        
        for (j, val) in values.iter().skip(1).enumerate() {
            if j >= num_kmers {
                break;
            }
            
            if val == &"1" {
                bitvec.set(j, true);
                hit_count += 1;
            }
        }
        
        if hit_count == 0 {
            zero_count += 1;
        }
        
        gene_hits.push(hit_count);
        matrix.push(bitvec);
    }
    
    // 输出统计信息
    println!("数据加载完成，耗时: {:?}", start.elapsed());
    println!("总基因数: {}", genes.len());
    println!("没有被任何k-mer覆盖的基因数: {}", zero_count);
    
    // 记录没有被覆盖的基因
    if zero_count > 0 {
        let uncovered_file = format!("{}_uncovered_genes.txt", output_prefix);
        let mut uncovered_writer = File::create(&uncovered_file)?;
        for (i, hits) in gene_hits.iter().enumerate() {
            if *hits == 0 {
                writeln!(uncovered_writer, "{}", genes[i])?;
            }
        }
        println!("未覆盖基因列表已保存到: {}", uncovered_file);
    }
    
    Ok(Data { genes, kmers, matrix, gene_hits })
}

fn load_sparse_matrix(file_path: &str, logger: &mut File, _output_prefix: &str) -> io::Result<Data> {
    let start = Instant::now();
    println!("开始读取稀疏矩阵文件...");
    writeln!(logger, "开始读取稀疏矩阵文件: {}", file_path)?;
    
    // 读取文件
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    
    // 从稀疏格式解析数据结构
    let mut gene_map: HashMap<String, usize> = HashMap::new();
    let mut kmer_map: HashMap<String, usize> = HashMap::new();
    let mut genes = Vec::new();
    let mut kmers = Vec::new();
    let mut sparse_data: HashMap<usize, HashSet<usize>> = HashMap::new();
    
    for (i, line_result) in reader.lines().enumerate() {
        let line = line_result?;
        
        // 跳过注释和空行
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 3 {
            writeln!(logger, "警告: 第{}行格式不正确: {}", i+1, line)?;
            continue;
        }
        
        let gene_id = parts[0].trim().to_string();
        let kmer_id = parts[1].trim().to_string();
        let value = parts[2].trim();
        
        // 只处理值为1的记录
        if value != "1" {
            continue;
        }
        
        // 添加基因和k-mer到映射
        let gene_idx = if let Some(&idx) = gene_map.get(&gene_id) {
            idx
        } else {
            let idx = genes.len();
            genes.push(gene_id.clone());
            gene_map.insert(gene_id, idx);
            idx
        };
        
        let kmer_idx = if let Some(&idx) = kmer_map.get(&kmer_id) {
            idx
        } else {
            let idx = kmers.len();
            kmers.push(kmer_id.clone());
            kmer_map.insert(kmer_id, idx);
            idx
        };
        
        // 添加到稀疏数据
        sparse_data.entry(gene_idx).or_insert_with(HashSet::new).insert(kmer_idx);
    }
    
    // 创建BitVec矩阵
    let mut matrix = Vec::with_capacity(genes.len());
    let mut gene_hits = Vec::with_capacity(genes.len());
    
    for gene_idx in 0..genes.len() {
        let mut bitvec = BitVec::from_elem(kmers.len(), false);
        let empty_set = HashSet::new();
        let kmer_indices = sparse_data.get(&gene_idx).unwrap_or(&empty_set);
        
        for &kmer_idx in kmer_indices {
            if kmer_idx < kmers.len() {
                bitvec.set(kmer_idx, true);
            }
        }
        
        gene_hits.push(kmer_indices.len());
        matrix.push(bitvec);
    }
    
    // 输出统计信息
    let zero_count = gene_hits.iter().filter(|&&count| count == 0).count();
    
    println!("数据加载完成，耗时: {:?}", start.elapsed());
    println!("总基因数: {}", genes.len());
    println!("总k-mer数: {}", kmers.len());
    println!("没有被任何k-mer覆盖的基因数: {}", zero_count);
    
    Ok(Data { genes, kmers, matrix, gene_hits })
}

fn greedy_kmer_selection(
    data: &Data, 
    target_counts: &[usize], 
    logger: &mut File, 
    _output_prefix: &str
) -> io::Result<Vec<(usize, usize, Vec<usize>)>> {
    let start = Instant::now();
    println!("开始贪心选择算法...");
    writeln!(logger, "开始贪心选择算法...")?;
    
    let num_genes = data.matrix.len();
    let num_kmers = data.kmers.len();
    let mut results = Vec::new();
    
    // 跟踪已覆盖的基因
    let mut covered_genes = BitVec::from_elem(num_genes, false);
    let mut selected = Vec::new();
    
    // 预先计算所有可能需要的k-mer
    let max_target = *target_counts.iter().max().unwrap_or(&1000);
    
    // 预计算k-mer覆盖信息
    println!("预计算k-mer覆盖信息...");
    let kmer_coverage: Vec<BitVec> = (0..num_kmers)
        .into_par_iter()
        .map(|kmer_idx| {
            let mut coverage = BitVec::from_elem(num_genes, false);
            for gene_idx in 0..num_genes {
                if data.matrix[gene_idx].get(kmer_idx).unwrap_or(false) {
                    coverage.set(gene_idx, true);
                }
            }
            coverage
        })
        .collect();
    
    // 添加进度条
    let pb = ProgressBar::new(max_target as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos}/{len} k-mers ({percent}%) {msg}")
        .unwrap()
        .progress_chars("=>-"));
    
    let mut current_covered = 0;
    
    while selected.len() < max_target && current_covered < num_genes {
        // 更新进度条
        pb.set_position(selected.len() as u64);
        let coverage_percent = (current_covered as f64 / num_genes as f64) * 100.0;
        pb.set_message(format!("覆盖率: {:.2}%", coverage_percent));
        
        // 并行计算每个候选k-mer的增益
        let best_kmer = (0..num_kmers)
            .into_par_iter()
            .filter(|&kmer_idx| !selected.contains(&kmer_idx))
            .map(|kmer_idx| {
                // 计算未覆盖的基因中，当前k-mer能覆盖的数量
                let mut gain = 0;
                for gene_idx in 0..num_genes {
                    if !covered_genes[gene_idx] && kmer_coverage[kmer_idx][gene_idx] {
                        gain += 1;
                    }
                }
                (kmer_idx, gain)
            })
            .reduce(|| (0, 0), |a, b| if a.1 > b.1 { a } else { b });
        
        if best_kmer.1 > 0 {
            let kmer_idx = best_kmer.0;
            selected.push(kmer_idx);
            
            // 更新已覆盖的基因
            for gene_idx in 0..num_genes {
                if !covered_genes[gene_idx] && kmer_coverage[kmer_idx][gene_idx] {
                    covered_genes.set(gene_idx, true);
                    current_covered += 1;
                }
            }
            
            // 检查是否达到目标数量
            for &target in target_counts {
                if selected.len() == target {
                    results.push((target, current_covered, selected.clone()));
                }
            }
        } else {
            break;
        }
    }
    
    pb.finish_with_message("k-mer选择完成");
    
    // 确保所有目标都有结果
    for &target in target_counts {
        if !results.iter().any(|(t, _, _)| *t == target) && !selected.is_empty() {
            let result_kmers = if selected.len() >= target {
                selected[0..target].to_vec()
            } else {
                selected.clone()
            };
            results.push((target, current_covered, result_kmers));
        }
    }
    
    // 按照target排序结果
    results.sort_by_key(|(target, _, _)| *target);
    
    // 输出最终结果
    let final_percent = (current_covered as f64 / num_genes as f64) * 100.0;
    println!("贪心算法完成，最终覆盖率: {}/{} ({:.2}%)，耗时: {:?}", 
             current_covered, num_genes, final_percent, start.elapsed());
    
    Ok(results)
} 