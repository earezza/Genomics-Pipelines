#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-jdilwort
#SBATCH --mem=32G
#SBATCH --array=0-2
#SBATCH --job-name=MEME-ChIP
#SBATCH --output=%x-%j.out

module load python/3.9 scipy-stack/2021a
source ngsENV/bin/activate

# Meme Suite
module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3
module load meme/5.4.1
module load bedtools/2.30.0


files=(SAMPLE_FILE_1_summits.bed SAMPLE_FILE_2_summits.bed SAMPLE_FILE_3_summits.bed)

# Extend summit range for peaks
bedtools slop -i Summits/${files[SLURM_ARRAY_TASK_ID]} -g Reference_Files/mm10/mm10.chrom.sizes -b 500 > Summits/${files[SLURM_ARRAY_TASK_ID]}_slop500.bed
sleep 60

# Get fasta sequence for extended summits
bedtools getfasta -fi Reference_Files/mm10/mm10.fa.masked -bed Summits/${files[SLURM_ARRAY_TASK_ID]}_slop500.bed -fo Summits/${files[SLURM_ARRAY_TASK_ID]}_slop500.bed.fasta
sleep 60

# Run meme-chip
meme-chip Summits/${files[SLURM_ARRAY_TASK_ID]}_slop500.bed.fasta -oc MEME-ChIP-OUTPUT_${files[SLURM_ARRAY_TASK_ID]}/ -ccut 500 -meme-nmotifs 15 -meme-mod anr -minw 15 -maxw 25 -db Reference_Files/HOCOMOCOv11_full_MOUSE_mono_meme_format.meme
