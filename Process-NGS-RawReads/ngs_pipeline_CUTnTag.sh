#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-jdilwort
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1G
#SBATCH --array=0-2
#SBATCH --job-name=CUTnTag_Process
#SBATCH --output=%x-%j.out

module load python/3.9 scipy-stack/2021a
module load samtools>=1.11 r>=4.0.5 bowtie2>=2.4.1 fastqc>=0.11.9
module load hisat2
module load java
source ngsENV/bin/activate
rm *=*

files=(SAMPLE_FOLDER_1 SAMPLE_FOLDER_2 SAMPLE_FOLDER_3)

python ngs_processing_pipeline.py -logfile ${files[SLURM_ARRAY_TASK_ID]}.log -reads ~/scratch/${files[SLURM_ARRAY_TASK_ID]}/ --technique cnt --species Mus --length 100 --reads_type paired --genome_index Reference_Files/bowtie2/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome --no_spikein -adapters 1 -outdir ~/scratch/${files[$SLURM_ARRAY_TASK_ID]}/
