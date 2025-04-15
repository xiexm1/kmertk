![kmertklogo](https://github.com/user-attachments/assets/5b18188a-79da-4001-ba04-7fb35b6f051e)
# KmerTK - K-mer工具集
[English](README.en.md) | 中文
KmerTK是一个高效的k-mer分析工具集，集成了多种k-mer相关的功能，专为高通量序列分析设计。

## 特性

- **Match**: 在序列中匹配k-mer并提取对应序列
- **Coverage**: 计算k-mer对序列的覆盖率
- **Distance**: 计算序列之间的编辑距离
- **Collection**: 构建k-mer矩阵和处理比对结果
- **MiniSet**: 选择最小k-mer子集来覆盖目标序列集合

## 安装

```bash
# 克隆仓库
git clone https://github.com/xiexm1/kmertk.git
cd kmertk

# 编译项目
cargo build --release
```

## 使用方法

### Match - k-mer匹配

```bash
kmertk match --kmer-file <kmer文件> --input-file <序列文件> --output-file <输出文件> [--extract-sequences]
```

### Coverage - 覆盖率分析

```bash
kmertk coverage --input-file <序列文件> --kmer-file <kmer文件> --output-file <输出文件> [--threshold <阈值>]
```

### Distance - 序列距离计算

```bash
kmertk distance --input-file <序列文件> --output-file <输出文件> [--pairs <计算对数>] [--all-pairs]
```

### Collection - k-mer矩阵构建

基于序列文件构建k-mer矩阵:
```bash
kmertk collection --input-mode seq --sequence-file <序列文件> --kmer-file <kmer文件> --output-file <输出文件> [--sparse] [--max-mismatches <错配数>]
```

基于SAM比对文件:
```bash
kmertk collection --input-mode sam --alignment-file <SAM文件> --output-file <输出文件> [--sparse] [--max-mismatches <错配数>]
```

基于Bowtie比对文件:
```bash
kmertk collection --input-mode bowtie --alignment-file <Bowtie文件> --output-file <输出文件> [--sparse]
```

### MiniSet - 最小k-mer集选择

```bash
kmertk miniset --matrix-file <矩阵文件> --output-prefix <输出前缀> [--sparse] [--threads <线程数>]
```

## 贡献

欢迎提交问题和拉取请求，一起改进KmerTK！如果觉得这个项目对你有帮助，欢迎给项目点个Star以示支持！

## 许可

本项目采用MIT许可证。 
