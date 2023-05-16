# xRead: an incremental overlapping graph construction approach
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/abPOA/blob/master/LICENSE)

## Table of Contents
1. [Getting started](#getting_stated)
2. [Introduction](#introduction)
3. [Installation](#installation)
    - [Building xRead from source files](#build_source)
    - [Downloading pre-built binary executable file](#build_binary)
4. [General usage](#general_usage)
    - [Input](#input)
    - [Output](#output)
    - [User cases](#user_cases)
    - [Commands and options](#commands_and_options)
5. [Contact](#contact)

## <a name="getting_stated"></a>Getting started

Downloading the source files from latest release:
```
wget https://github.com/tcKong47/xRead/releases/download/v1.0.0/xRead-v1.0.0.tar.gz
tar -zxvf xRead-v1.0.0.tar.gz && cd xRead-v1.0.0
```
Building from source files and running with the test data:
```
make
./xRead ./test_data/reads.fasta -f test_data
```

## <a name="introduction"></a>Introduction
xRead is an incremental overlapping graph construction approach that is able to achieve high scalability, performance, and yields simultaneously. Guided by a novel read-depth-based objective function, xRead iteratively builds and refines the overlapping graph with heuristic read indexing and lightweight alignment skeletons and output the core overlap information in PAF format.

xRead has outstanding scalability for memory usage which enables to build the overlapping graphs of high-coverage sequencing datasets of very large genomes with low and tunable RAM space cost, e.g., building an overlapping graph of 32X PacBio sequencing dataset (1.9 Terabyte) of Axolotl genome with 64GB RAM. Moreover, xRead has high speed for various genome sizes from E. coli to human genomes, and can enable the production of highly accurate read overlaps without loss of sensitivity, which helps to construct high-quality overlapping graphs. xRead is suited to handle large genomes/datasets produced by ONT and PacBio platforms with only small or medium servers/clusters and be cost-effective large-scale de novo assembly tasks having numbers of genomes to be assembled. 

## <a name="installation"></a>Installation
### <a name="build_source"></a>Building xRead from source files
You can build xRead from source files. The source code of xRead is written in C, and has been tested on 64-bit Linux. The source code can be downloaded from the latest release or git clone command. Please make sure you have gcc (>=6.4.0) installed before compiling.
```
wget https://github.com/tcKong47/xRead/releases/download/v1.4.1/xRead-v1.0.0.tar.gz
tar -zxvf xRead-v1.0.0.tar.gz
cd xRead-v1.0.0; make
```
```
git clone --recursive https://github.com/tcKong47/xRead.git
cd xRead; make
```

### <a name="build_binary"></a>Downloading pre-built binary executable file
Or you can download the pre-built binary file if you meet compiling issues.
```
wget https://github.com/tcKong47/xRead/releases/download/v1.0.0/xRead-v1.0.0_x64-linux.tar.gz
tar -zxvf xRead-v1.0.0_x64-linux.tar.gz
```

## <a name="general_usage"></a>General usage
### <a name="input"></a>Input
xRead works with FASTA, FASTQ, gzip'd FASTA (.fa.gz), and gzip'd FASTQ (.fq.gz) formats. The input file is expected to contains all sequences of a single sample of a specific organism to perform the incremental overlap graph construction.

xRead also works with PAF and gzip'd (.paf.gz) formats and provides an additional function to expand the initial core graphs to comprehensive graphs based on the further analysis of transitive relationships among the overlaps. If the option *-N --trans-file* is enabled, xRead skips the overlapping process and takes the PAF file as input to directly performed the expanding iterations, which the number of iterations is defined by option *-R --trans-iter*.

### <a name="output"></a>Output
xRead generates a collection of overlaps between input reads/overlap relationships and outputs them in the PAF format. For more information about PAF format, please refer to PAF: a Pairwise mApping Format (https://github.com/lh3/miniasm/blob/master/PAF.md).

### <a name="user_cases"></a>User cases
To construct overlapping graph and output in PAF format.
```
xRead sequence.fa/fq -f output_path/output_filename_prefix
```

To construct overlapping graph and expand it to comprehensive overalpping graph via transitive searching of -R iterations.
```
xRead sequence.fa/fq -f output_path/output_filename_prefix -R 3
```

To expand overlapping graph to comprehensive overalpping graph.
```
xRead overlaps.paf -NR 3
```

### <a name="commands_and_options"></a>Commands and options
***Usage of xRead***
```
Usage: xRead <reads.fa/fq>/<overlaps.paf> [options]

    <reads.fq/fa>           The input reads in fasta or fastq format, necessary if -N is not enabled by users.
    <overlaps.paf>          The input overlapping graph in paf format, necessary if -N is enabled to skip the overlapping discovery step.
```

***Program options***
| Short option | Long option | Description | Default |
| :----------: | :---------- | :---------- | :-----------: |
| -h           | --help | Print help menu. | NULL |
| -f           | --output-file | The path and name prefix of output file. | NULL |
| -M           | --memory | Maximum allowed memory. | 16 |
| -t           | --thread-n | Number of used threads. | 8 |
| -p           | --read-type | Specify the type of reads to calculating the warting length for alignment <br> skeleton.<br> *1*: For reads with high error rates ~15%.<br> *2*: For reads with low error rates ~1%. | 1 |

***Algorithm options***
| Short option | Long option | Description | Default |
| :----------: | :---------- | :---------- | :-----------: |
| -x | --x-longest-read | The longest x% reads indexed for the first iteration. | 3.0 |
| -X | --x-random-read | The longest X% reads indexed for the second and later iteration. | 10.0 |
| -k | --k-mer | The length of k-mer for sequences sketching. | 15 |
| -w | --window-size | Window size for sequences sketching. | 5 |
| -l | --l-mer | The length of l-mer of the auxiliary index hash table. | 11 |
| -r | --repeat-n | The proportion for filtering the most repetitive minimizer hits. | 0.005 |
| -e | --search-step | The number of search step of SDP graph construction. | 10 |
| -n | --top-n | The number of overlaps each read retained for each iteration. | 2 |
| -b | --matching-bases | Minimum required matching bases for a single overlap between two reads. | 100 |
| -m | --min-ove | Minimum required length of overlap. | 500 |
| -a | --max-hang | Maximum allowed total length of left and right overhangs. | 2000 |
| -L | --split-len | The length of split blocks for read coverage estimation. | 1000 |
| -S | --read-part-size |  Estimating the average coverage of every S consecutive blocks to mark the <br> less-covered read portions. | 3 |
| -c | --cov-ratio | The proportion of less-covered reads for next iteration. | 0.5 |
| -I | --iter-times | Maximum allowed number of iterations. | 8 |
| -R | --trans-iter | Enabling the calculating of transitive overlaps for R iterations. | 0 |
| -N | --trans-only | Skipping the overlapping discovery, only perfomred R transitive iterations to <br> generate comprehensive graph. | NULL |

***Recommended options***

For ONT and PacBio CLR data:
```
-k/--k-mer       15
-w/--window-size  5
```
For PacBio CCS (HiFi) and ONT data of super high accuracy mode:
```
-k/--k-mer       19
-w/--window-size 40
```

## <a name="contact"></a>Contact
Please post on [github issues](https://github.com/tcKong47/xRead/issue) or contact Tangchao Kong 21B903020@stu.hit.edu.cn for questioning, advising, and reporting bugs.
