## Description  
Pipeline to process raw reads from next-generation sequencing experiments (CUT&Tag, ChIP-Seq, ATAC-Seq, RNA-Seq, MNase-Seq, etc...).  

## Usage
Place .fastq.gz files for a sample replicate into a folder.  

&emsp;e.g. paired-end reads:  
&emsp;(Note: one paired-end replicate will have an "_R1" and "_R2" file)  
&emsp;&emsp;RAW_READS/  
&emsp;&emsp;&emsp;read_R1.fastq.gz  
&emsp;&emsp;&emsp;read_R2.fastq.gz  

&emsp;Run script (see default options):  
&emsp;e.g. RNA-Seq:  
> python ngs_processing_pipeline.py -logfile RAW_READS.log -reads RAW_READS/ --technique rnaseq --species Mus --length 100 --reads_type paired --genome_index hisat_genomes_index/UCSC/mm10/genome --no_spikein -adapters 1 -qctrim -outdir RAW_READS/  

&emsp;e.g. Other:  
> python ngs_processing_pipeline.py -logfile RAW_READS.log -reads RAW_READS/ --species Mus --length 100 --reads_type paired --genome_index mm10/Bowtie2Index/genome --no_spikein -adapters 1 -qctrim -outdir RAW_READS/  

&emsp;e.g. single-end reads:  
&emsp;&emsp;RAW_READS/  
&emsp;&emsp;&emsp;read_1.fastq.gz  
&emsp;&emsp;&emsp;read_2.fastq.gz  
&emsp;&emsp;&emsp;read_3.fastq.gz  
&emsp;&emsp;&emsp;read_4.fastq.gz  

## Software Required  
FastQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
MultiQC https://multiqc.info/  
Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  
Hisat2 https://daehwankimlab.github.io/hisat2/  
Picard https://broadinstitute.github.io/picard/  
Samtools  http://www.htslib.org/  
Deeptools  https://deeptools.readthedocs.io/en/develop/  
SEACR https://github.com/FredHutch/SEACR  
MACS https://github.com/macs3-project/MACS    
GoPeaks https://github.com/maxsonBraunLab/gopeaks  
  
## Reference Genome Files  
Bowtie2 <a href="https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">genome index files</a>  for alignment (on right under Indexes). Alternatively, you can download from <a href="https://support.illumina.com/sequencing/sequencing_software/igenome.html">iGenomes.</a>  
Hisat2 <a href="https://daehwankimlab.github.io/hisat2/download/">genome index files</a> for RNA-Seq alignment.  

### Blacklist Regions  
If blacklisted regions wish to be removed in bamCoverage, you can find these files <a href="https://github.com/Boyle-Lab/Blacklist">here.</a>  
___  

### Additional Resources  
CUT&Tag https://www.nature.com/articles/s41467-019-09982-5  
Pipeline example https://yezhengstat.github.io/CUTTag_tutorial/  
Additional https://learn.gencore.bio.nyu.edu/  
  
<img src="ngs_pipeline.png" alt="Pipeline Process">    

___     
