use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use flate2::read::GzDecoder;

/// 解析FASTA文件，支持常规文件和gzip压缩文件
/// 返回 Vec<(序列ID, 序列)>
pub fn parse_fasta<P: AsRef<Path>>(file_path: P) -> io::Result<Vec<(String, String)>> {
    let path = file_path.as_ref();
    let reader: Box<dyn io::Read> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(path)?))
    } else {
        Box::new(File::open(path)?)
    };
    
    let buf_reader = BufReader::new(reader);
    let mut sequences = Vec::new();
    let mut current_id = String::new();
    let mut current_seq = String::new();

    for line in buf_reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() {
            continue;
        } else if line.starts_with('>') {
            if !current_id.is_empty() {
                sequences.push((current_id.clone(), current_seq.clone()));
                current_seq = String::new();
            }
            // 提取序列ID，只保留'>'后的第一个单词
            current_id = line[1..].split_whitespace().next()
                .unwrap_or("").to_string();
        } else {
            current_seq.push_str(line);
        }
    }
    
    if !current_id.is_empty() {
        sequences.push((current_id, current_seq));
    }
    
    Ok(sequences)
}

/// 将FASTA序列转换为HashMap格式，方便通过ID快速查找
pub fn fasta_to_map<P: AsRef<Path>>(file_path: P) -> io::Result<HashMap<String, String>> {
    let sequences = parse_fasta(file_path)?;
    Ok(sequences.into_iter().collect())
}

/// 读取FASTA文件中的序列数量
pub fn count_sequences<P: AsRef<Path>>(file_path: P) -> io::Result<usize> {
    let path = file_path.as_ref();
    let reader: Box<dyn io::Read> = if path.to_string_lossy().ends_with(".gz") {
        Box::new(GzDecoder::new(File::open(path)?))
    } else {
        Box::new(File::open(path)?)
    };
    
    let buf_reader = BufReader::new(reader);
    let mut count = 0;

    for line in buf_reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            count += 1;
        }
    }
    
    Ok(count)
} 