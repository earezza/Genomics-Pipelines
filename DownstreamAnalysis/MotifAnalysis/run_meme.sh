#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --account=def-jdilwort
#SBATCH --mem=32G
#SBATCH --job-name=MEME
#SBATCH --output=%x-%j.out

module load python/3.9 scipy-stack/2021a
source ngsENV/bin/activate

# Meme Suite
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load meme/5.4.1
module load bedtools/2.30.0


# PREPARE QUERY SEQUENCES
# Extend summit range for peaks
bedtools slop -i Summits/query_summits.bed -g Reference_Files/mm10/mm10.chrom.sizes -b 500 > Summits/query_summits_slop500.bed
# Get fasta sequence for extended summits
bedtools getfasta -fi Reference_Files/mm10/mm10.fa.masked -bed Summits/query_summits_slop500.bed -fo Summits/query_summits_slop500.bed.fasta


# PREPARE CONTROL SEQUENCES
# Extend summit range for peaks
bedtools slop -i Summits/Gene_list_sorted.bed -g Reference_Files/mm10/mm10.chrom.sizes -b 500 > Summits/Gene_list_sorted_slop500.bed
# Get fasta sequence for extended summits
bedtools getfasta -fi Reference_Files/mm10/mm10.fa.masked -bed Summits/Gene_list_sorted_slop500.bed -fo Summits/Gene_list_sorted_slop500.bed.fasta

# Run MEME
meme Summits/query_summits_slop500.bed.fasta -oc MEME-OUTPUT/ -objfun de -neg Summits/Gene_list_sorted_slop500.bed.fasta -mod anr -minw 15 -maxw 25 -nmotifs 15 -dna -revcomp
