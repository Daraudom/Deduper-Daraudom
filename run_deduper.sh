#!/usr/bin/env bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --job-name=deduper.sh
#SBATCH --output=deduper_output.log
#SBATCH --error=deduper_error.log

# MUST TAKE A SORTED SAM FILE
# uncomment below when needed
/usr/bin/time -v ./deduper.py -f input/sorted_test_big.sam -o big_test.sam -u STL96.txt 
#./deduper.py -f Unit_Tests/input_softclip_1.sam -o Unit_Tests/output_softclip_1.sam -u STL96.txt 

