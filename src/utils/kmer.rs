use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// 从文件中读取k-mer
pub fn read_kmers<P: AsRef<Path>>(kmer_file: P) -> io::Result<Vec<String>> {
    let file = File::open(kmer_file)?;
    let reader = BufReader::new(file);
    let kmers: Vec<String> = reader
        .lines()
        .filter_map(|line| line.ok())
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty())
        .collect();
    
    Ok(kmers)
}

/// 查找字符串中所有k-mer的位置
pub fn find_kmer_positions(sequence: &str, kmer: &str) -> Vec<usize> {
    let mut positions = Vec::new();
    let mut start = 0;
    
    while let Some(pos) = sequence[start..].find(kmer) {
        let abs_pos = start + pos;
        positions.push(abs_pos);
        start = abs_pos + 1;
    }
    
    positions
}

/// 找出字符串中所有k-mer匹配的区间
pub fn find_kmer_intervals(sequence: &str, kmers: &HashSet<&str>) -> Vec<(usize, usize)> {
    if kmers.is_empty() {
        return Vec::new();
    }
    
    let k = kmers.iter().next().unwrap().len();
    
    kmers.iter()
        .flat_map(|kmer| {
            find_kmer_positions(sequence, kmer)
                .into_iter()
                .map(move |pos| (pos, pos + k - 1))
        })
        .collect()
}

/// 合并重叠的区间
pub fn merge_intervals(mut intervals: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if intervals.is_empty() {
        return intervals;
    }
    
    intervals.sort_by_key(|&(start, _)| start);
    let mut merged = Vec::new();
    let mut current = intervals[0];
    
    for &(start, end) in intervals.iter().skip(1) {
        if start <= current.1 + 1 {
            // 区间重叠或相邻，合并
            current.1 = current.1.max(end);
        } else {
            // 区间不重叠，保存当前区间并开始新区间
            merged.push(current);
            current = (start, end);
        }
    }
    
    merged.push(current);
    merged
}

/// 计算序列被k-mers覆盖的比例
pub fn calculate_coverage(sequence: &str, kmers: &HashSet<&str>) -> f64 {
    let intervals = find_kmer_intervals(sequence, kmers);
    let merged = merge_intervals(intervals);
    
    let covered_length: usize = merged.iter()
        .map(|&(start, end)| end - start + 1)
        .sum();
    
    let seq_len = sequence.len();
    if seq_len == 0 {
        return 0.0;
    }
    
    covered_length as f64 / seq_len as f64
}

/// 计算两个序列之间的Levenshtein距离
pub fn levenshtein(a: &str, b: &str) -> usize {
    let a_len = a.chars().count();
    let b_len = b.chars().count();
    
    // 边界情况
    if a_len == 0 { return b_len; }
    if b_len == 0 { return a_len; }
    
    let mut dp = vec![vec![0; b_len + 1]; a_len + 1];
    
    // 初始化第一行和第一列
    for i in 0..=a_len {
        dp[i][0] = i;
    }
    for j in 0..=b_len {
        dp[0][j] = j;
    }
    
    // 计算编辑距离
    for (i, a_char) in a.chars().enumerate() {
        for (j, b_char) in b.chars().enumerate() {
            let cost = if a_char == b_char { 0 } else { 1 };
            dp[i + 1][j + 1] = std::cmp::min(
                dp[i][j + 1] + 1,           // 删除
                std::cmp::min(
                    dp[i + 1][j] + 1,       // 插入
                    dp[i][j] + cost         // 替换或匹配
                )
            );
        }
    }
    
    dp[a_len][b_len]
} 