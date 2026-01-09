# Deduper

A memory-efficient deduplication and Python-based tool for removing PCR duplicates from single-end SAM files using known UMIs.

## Overview

PCR amplification during sequencing library preparation can introduce duplicate reads that artificially inflate signal. This tool removes PCR duplicates from **single-end, uniquely mapped SAM files** by identifying reads that share the same alignment position, strand orientation, and UMI.

## Usage

```bash
python deduper.py \
  --file sorted_input.sam \
  --umi known_umis.txt \
  --output deduplicated.sam
```
## Inputs

- `-f --file`: A **single-end sorted `.sam` file** containing uniquely mapped reads
- `-u --umi`: A text file of known UMIs (Unique Molecular Identifiers)
- `-o --output`: A name of the output file to hold the deduplicated reads.


## What Is Considered a Duplicate?

Two reads are considered PCR duplicates if **all** of the following conditions are met:

### 1. Same Alignment Position
- **Same chromosome**  
  → SAM column 3 (`RNAME`)
- **Same strand**  
  → Determined from SAM column 2 (`FLAG`)
  - `FLAG & 16 != 0` → reverse strand (`-`)
  - otherwise → forward strand (`+`)
- **Same true 5′ start position**
  - Derived from:
    - SAM column 4 (`POS`)
    - SAM column 6 (`CIGAR`)
  - Soft clipping (`S`) in the CIGAR string is used to adjust the 5′ position appropriately depending on strand orientation

### 2. Same UMI
- The UMI must:
  - Match exactly
  - Be present in the provided UMI whitelist

## Read Filtering Rules

The following reads are excluded from processing:

- Reads with **unknown UMIs** (not found in the whitelist)
- **Unmapped reads** (`FLAG & 4 != 0`)
- PCR duplicates beyond the first observed instance

