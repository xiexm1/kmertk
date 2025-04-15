use clap::{Parser, Subcommand};
use std::io;

mod utils;
mod modules;

use modules::{
    match_mod,
    distance,
    coverage,
    miniset,
    collection,
};

/// KmerTK: 一体化k-mer分析工具包
#[derive(Parser)]
#[command(name = "kmertk")]
#[command(author = "KmerTK Team")]
#[command(version = "0.1.0")]
#[command(about = "K-mer分析工具集", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// 匹配k-mer并提取匹配的序列
    Match(match_mod::MatchArgs),
    
    /// 计算序列间的距离
    Distance(distance::DistanceArgs),
    
    /// 计算k-mer的覆盖率
    Coverage(coverage::CoverageArgs),
    
    /// 选择最小的k-mer集合覆盖基因
    MiniSet(miniset::MiniSetArgs),
    
    /// 构建k-mer矩阵和处理序列比对
    Collection(collection::CollectionArgs),
}

fn main() -> io::Result<()> {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Match(args) => match_mod::run_match(&args),
        Commands::Distance(args) => distance::run_distance(&args),
        Commands::Coverage(args) => coverage::run_coverage(&args),
        Commands::MiniSet(args) => miniset::run_miniset(&args),
        Commands::Collection(args) => collection::run_collection(&args),
    }
}
