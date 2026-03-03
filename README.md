# Deduper 🧬  
UMI-Aware PCR Duplicate Removal for Sorted SAM Files (Single-End)
Deduper is a Python command-line tool that removes PCR duplicates from a coordinate-sorted SAM file using Unique Molecular Identifiers (UMIs).

A read is considered unique based on the following profile:
(UMI, Chromosome, Strand, Corrected 5′ Position)
Only the first occurrence of each unique profile is retained. All subsequent duplicates are discarded.


## Overview

PCR amplification can introduce duplicate reads that inflate downstream analyses. Deduper removes these duplicates using UMI-aware logic and strand-aware 5′ coordinate correction derived from the CIGAR string.

The program:
- Reads a sorted SAM file
- Extracts UMIs from read names
- Validates UMIs against a whitelist
- Determines strand using the SAM FLAG field
- Corrects the 5′ mapping position using CIGAR operations
- Writes a deduplicated SAM output file
- Generates a summary statistics report (`report.txt`)

## Requirements

- Python 3.9+
- NumPy

Install NumPy if needed:

conda install numpy


## Input Requirements

### 1. Sorted SAM File

The input SAM file must be sorted by chromosome and position.

Example:
```
samtools sort -O SAM input.sam -o input.sorted.sam
```

### 2. Single-End Reads

This script is designed for single-end sequencing data.

### 3. UMI Format

UMIs must appear at the end of the QNAME field, separated by a colon `:`.
```
Example:

READ123:ATCG

Extracted UMI → ATCG
```

### 4. UMI Whitelist File

Provide a text file containing valid UMIs (one per line).

Example:

ATCG
GGTA
CCTA
...

Reads containing UMIs not present in this file are ignored.


## Usage
```
./deduper.py -f input.sorted.sam -u umi_list.txt -o deduped.sam
```

### Arguments

- `-f / --file`   (required)  Input sorted SAM file
- `-u / --umi`    (required)  File containing valid UMIs
- `-o / --output` (optional)  Output SAM filename (default: deduper_output)


## How Duplicate Detection Works

For each non-header SAM record:

1. Parse:
   - QNAME
   - FLAG
   - RNAME (chromosome)
   - POS
   - CIGAR

2. Extract:
   - UMI from QNAME (substring after last `:`)
   - Strand from FLAG  
     + if bit 16 not set  
     - if bit 16 is set  

3. Compute corrected 5′ position:
   - On + strand → adjust for leading soft clipping
   - On - strand → compute 5′ using alignment span derived from CIGAR

4. Construct unique key:

(UMI, chromosome, strand, corrected_5prime_position)

5. If key already exists → duplicate (discard)  
   If new → write read to output and store key


## Output Files

### Deduplicated SAM File
- All header lines are preserved
- Only unique reads are retained

### report.txt

Contains:
- Total number of reads processed
- Total duplicated reads removed
- Total unique reads retained
- Per-chromosome unique read counts

Example:

Deduper Report
Total Number of Reads: 500000
Total Number of Duplicated Reads: 120000
Total Number of Unique Reads: 380000
chr1: 90000
chr2: 85000


## Limitations

- Requires coordinate-sorted SAM input
- Designed for single-end reads
- Assumes UMI appears at end of QNAME
- Does not support paired-end deduplication
- Memory usage depends on number of unique read profiles per chromosome
