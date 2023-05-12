# xRead: an incremental overlapping graph construction approach
[![License](https://img.shields.io/badge/License-MIT-black.svg)](https://github.com/yangao07/abPOA/blob/master/LICENSE)

## Table of Contents
1. [Getting started](#getting_stated)
2. [Introduction](#introduction)
3. [Installation](#installation)
    - Building xRead from source files
    - Downloading Pre-built binary executable file
4. [General usage](#general_usage)
    - Input
    - Output
    - Commands and options
    - User cases

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
xRead is an incremental overlapping graph construction approach that is able to achieve high scalability, performance, and yields simultaneously. Guided by a novel read-depth-based objective function, xRead iteratively builds and refines the overlapping graph with heuristic read indexing and lightweight alignment skeletons. 

xRead has outstanding scalability for memory usage which enables to build the overlapping graphs of high-coverage sequencing datasets of very large genomes with low and tunable RAM space cost, e.g., building an overlapping graph of 32X PacBio sequencing dataset (1.9 Terabyte) of Axolotl genome with 64GB RAM. Moreover, xRead has high speed for various genome sizes from E. coli to human genomes, and can enable the production of highly accurate read overlaps without loss of sensitivity, which helps to construct high-quality overlapping graphs. 

xRead is suited to handle large genomes/datasets produced by ONT and PacBio platforms with only small or medium servers/clusters and be cost-effective large-scale de novo assembly tasks having numbers of genomes to be assembled. 

For more information, please refer to our paper published in xxx.

## <a name="installation"></a>Installation
### Building xRead from source files
You can build xRead from source files. The source code of xRead is written in C, and has been tested on 64-bit Linux. The source code can be downloaded from the latest release or git clone command. Please make sure you have gcc (>=6.4.0) installed before compiling.
```
wget https://github.com/tcKong47/xRead/releases/download/v1.4.1/xRead-v1.0.0.tar.gz
tar -zxvf xRead-v1.0.0.tar.gz
cd xRead-v1.0.0; make

git clone --recursive https://github.com/tcKong47/xRead.git
cd xRead; make
```

### Downloading Pre-built binary executable file
Or you can download the pre-built binary file if you meet compiling issues.
```
wget https://github.com/tcKong47/xRead/releases/download/v1.0.0/xRead-v1.0.0_x64-linux.tar.gz
tar -zxvf xRead-v1.0.0_x64-linux.tar.gz
```

## <a name="general_usage"></a>General usage
### Input
xRead works with FASTA, FASTQ, gzip'd FASTA(.fa.gz) and gzip'd FASTQ(.fq.gz) formats. The input file is expected to contains all sequences of a single sample of a specific organism to perform the incremental overlapping graph construction.

### Output
xRead generates a collection of overlaps between input reads and outputs them in the PAF format. For more information about PAF format, please refer to PAF: a Pairwise mApping Format(https://github.com/lh3/miniasm/blob/master/PAF.md).

### Commands and options

### User cases
