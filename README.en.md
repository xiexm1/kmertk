# KmerTK - K-mer Toolkit

KmerTK is an efficient k-mer analysis toolkit that integrates multiple k-mer related functionalities, designed for high-throughput sequence analysis.

## Features

- **Match**: Match k-mers in sequences and extract corresponding sequences
- **Coverage**: Calculate k-mer coverage of sequences
- **Distance**: Calculate edit distance between sequences
- **Collection**: Build k-mer matrices and process alignment results
- **MiniSet**: Select minimal k-mer subsets to cover target sequence sets

## Installation

```bash
# Clone the repository
git clone https://github.com/xiexm1/kmertk.git
cd kmertk

# Build the project
cargo build --release
```

## Usage

### Match - K-mer Matching

```bash
kmertk match --kmer-file <kmer_file> --input-file <sequence_file> --output-file <output_file> [--extract-sequences]
```

### Coverage - Coverage Analysis

```bash
kmertk coverage --input-file <sequence_file> --kmer-file <kmer_file> --output-file <output_file> [--threshold <threshold>]
```

### Distance - Sequence Distance Calculation

```bash
kmertk distance --input-file <sequence_file> --output-file <output_file> [--pairs <number_of_pairs>] [--all-pairs]
```

### Collection - K-mer Matrix Construction

Build k-mer matrix based on sequence file:
```bash
kmertk collection --input-mode seq --sequence-file <sequence_file> --kmer-file <kmer_file> --output-file <output_file> [--sparse] [--max-mismatches <mismatches>]
```

Using SAM alignment file:
```bash
kmertk collection --input-mode sam --alignment-file <SAM_file> --output-file <output_file> [--sparse] [--max-mismatches <mismatches>]
```

Using Bowtie alignment file:
```bash
kmertk collection --input-mode bowtie --alignment-file <Bowtie_file> --output-file <output_file> [--sparse]
```

### MiniSet - Minimal K-mer Set Selection

```bash
kmertk miniset --matrix-file <matrix_file> --output-prefix <output_prefix> [--sparse] [--threads <number_of_threads>]
```

## Contributing

Issues and pull requests are welcome to help improve KmerTK! If this project is helpful to you, please consider giving it a star to show your support!

## License

This project is licensed under the MIT License.
