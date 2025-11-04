#!/usr/bin/env python

import argparse
import re
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description="Given a sorted sam file of single-end reads, remove \
                                     PCR duplicates.")
    parser.add_argument("-f", "--file", help="Provide a sorted sam file", required=True)
    parser.add_argument("-o", "--output", help="Specify an output file. Default is set to \
                        deduper_output.sam", default="deduper_output")
    parser.add_argument("-u", "--umi", help="Specify the input file that contains the list of UMIs",
                        required= True)
    
    return parser.parse_args()

# Set global variables
args = get_args()
input_file = args.file
output_file = args.output
umi_file = args.umi

# stats variables
num_reads = 0
duplicated_reads = 0

# -------DEFINE HELPER FUNCTIONS-----------

def parse_umi_file(umi: str):
    """
    Take the UMI txt file and return a string list of the 96 UMIs
    """
    umi_list = []
    with open(umi, 'r') as fin:
        for line in fin:
            umi_list.append(line.strip())
    return umi_list

def parse_sam_file(sam_record: str):
    """
    Input: A sam file record
    Output: A tuple that contains:
        - UMI
        - Chromosome Number
        - Left-Most Mapped Position
        - Strandness
        - CIGAR string
    """
    # Break the sam record into its columns
    sam_record_split = sam_record.split(sep='\t')

    query_col = sam_record_split[0]
    bitwise_flg = int(sam_record_split[1])
    chr = sam_record_split[2]
    left_pos = int(sam_record_split[3])
    cigar = sam_record_split[5]

    # extract umi
    umi = query_col.split(sep=':')[-1]

    # extract strandness
    strand = "+"
    if (bitwise_flg & 16) == 16:
        strand = "-"

    return (umi, chr, left_pos, strand, cigar)

def parse_cigar_str(cigar: str, strand: str):
    modified_cig = []

    cig_split = re.findall("[0-9]*[A-Z]", cigar)

    for i in range(len(cig_split)):
        # get the regex matches
        num = re.findall("[0-9]*", cig_split[i])[0]
        if i == 0:
            if "S" in cig_split[i]: # with soft clip
                if strand == "-": continue
                # add a minus sign for soft clip and +
                modified_cig.append(-int(num))
            elif strand == "+": # no soft clip and +
                modified_cig.append(0)
                break
            else: # no soft clip and -
                modified_cig.append(int(num))
        elif "I" in cig_split[i]:
            modified_cig.append(0)
        else: # If it's M or D or X or S in the end
            modified_cig.append(int(num))
    return modified_cig

def calc_five_prime(left_pos: int, mod_cigar: list[int], strand: str):
    true_pos = 0

    if strand == "+":
        true_pos = left_pos + mod_cigar[0]
    else:
        true_pos = left_pos + np.sum(mod_cigar)
    
    return true_pos

# Parse the UMI File
umi_list = parse_umi_file(umi_file)

# Define global set
profile_set = set()

# testing variables
chr_read_count = {}

# Parse the sam file
with open (input_file, 'r') as sam_lines, \
    open (output_file, 'w') as fout:
    # initialize current chromosome
    curr_chr = None
    
    while True:
        sam_line = sam_lines.readline()

        # EOF , reset profile and update
        if sam_line == '':
            chr_read_count[curr_chr] = len(profile_set)
            duplicated_reads += len(profile_set)

            # actual stuff
            profile_set.clear()
            curr_chr = chr

            break
        # if starts with @ just write it out to output file
        if sam_line.startswith("@"):
            fout.write(sam_line)
            continue

        umi, chr, left_pos, strand, cigar = parse_sam_file(sam_line)
        num_reads += 1 # increment read

        # base case set
        if curr_chr == None:
            curr_chr = chr
        # for transitioning
        if curr_chr != chr:
            chr_read_count[curr_chr] = len(profile_set)
            duplicated_reads += len(profile_set)
            # actual stuff
            profile_set.clear()
            curr_chr = chr
            
        # validate umi
        if umi not in umi_list: continue

            # checking chr1

        # calculate 5 pos
        mod_cig_line = parse_cigar_str(cigar=cigar, strand=strand)
        five_pos = calc_five_prime(left_pos=left_pos, \
                                   mod_cigar=mod_cig_line, \
                                    strand=strand)
        
        # create profile 
        profile_sam = (umi, chr, strand, five_pos)

        # debug
        #print(profile_sam)

        if profile_sam not in profile_set:
            profile_set.add(profile_sam)
            fout.write(sam_line)

        

#print(chr_read_count)
print(f"total num of unique reads: {sum(list(chr_read_count.values()))}")

with open("report.txt", "w") as stats_file:
    stats_file.write("Deduper Report\n")
    stats_file.write(f"Total Number of Reads: {num_reads}\n")
    stats_file.write(f"Total Number of Duplicated Reads: {duplicated_reads}\n")
    stats_file.write(f"Total Number of Unique Reads: {num_reads-duplicated_reads}\n")

    for chr, val in chr_read_count.items():
        stats_file.write(f'{chr}: {val}\n')
