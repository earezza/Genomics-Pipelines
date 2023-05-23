#!/bin/bash
# After logging in, enter -> cd projects/def-jdilwort/$USER/
# Then the script below will prep everything and may take several minutes
# Note that only mm10 and hg38 genome assemblies are downloaded automatically, others require manual downloading

date

# --- Download GitHub ---
git clone https://github.com/Dilworth-Laboratory/Genomics-Pipelines.git
cp Genomics-Pipelines/Process-NGS-RawReads/ngs_processing_pipeline.py .
cp Genomics-Pipelines/Process-NGS-RawReads/ngs_pipeline_*.sh .


# --- Download external software ---
wget -O picard.jar https://github.com/broadinstitute/picard/releases/download/3.0.0/picard.jar
git clone https://github.com/FredHutch/SEACR.git
wget -O gopeaks https://github.com/maxsonBraunLab/gopeaks/releases/download/v1.0.0/gopeaks-linux-amd64
chmod +x gopeaks

# --- Create virtualenvironment ---
module load python/3.9
virtualenv --no-download ngsENV
source ngsENV/bin/activate
pip install --no-index --upgrade pip
pip install -r Genomics-Pipelines/Process-NGS-RawReads/cc_requirements.txt


# --- Download reference genomes ---
mkdir Reference_Files/
cd Reference_Files/

mkdir bowtie2/
cd bowtie2/
#wget -O Homo_sapiens_UCSC_hg38.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
#tar -xzvf Homo_sapiens_UCSC_hg38.tar.gz
wget -O Mus_musculus_UCSC_mm10.tar.gz http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz
tar -xzvf Mus_musculus_UCSC_mm10.tar.gz
cd ../

mkdir hisat2/
cd hisat2/
#wget -O hg38_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz
#tar -xzvf hg38_genome.tar.gz
wget -O mm10_genome.tar.gz https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
tar -xzvf mm10_genome.tar.gz
cd ../../

date
