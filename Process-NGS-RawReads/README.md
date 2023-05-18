## Description  
Pipeline to process raw reads from next-generation sequencing experiments (CUT&Tag, ChIP-Seq, ATAC-Seq, RNA-Seq, MNase-Seq, etc...).  
<strong>This should be run on an HPC platform such as ComputeCanada.</strong>  
___  
## Initial Setup  
#### These steps need to be performed once:  
Note: Automatic setup if you only download (or copy/paste) and run ./ngs_setup.sh which takes several minutes to complete.
<ol>  
  <li>Download/clone this https://github.com/earezza/Genomics-Pipelines to your working directory</li>  
  <li>Load Python  
    <ul>  
      <li>module load python/3.9</li>  
    </ul>  
  </li>  
  <li>Create virtualenv, activate, and install packages from cc_requirements.txt  
    <ul>  
      <li>virtualenv --no-download ngsENV</li>  
    </ul>  
    <ul>  
      <li>source ngsENV/bin/activate</li>  
    </ul>  
    <ul>  
      <li>pip install --no-index --upgrade pip</li>  
    </ul>  
    <ul>  
      <li>pip install -r cc_requirements.txt</li>  
    </ul>  
  </li>  
  <li>Download Picard https://broadinstitute.github.io/picard/ to your working directory</li>  
  <li>Download/clone SEACR https://github.com/FredHutch/SEACR to your working directory</li>  
  <li>Download/clone GoPeaks https://github.com/maxsonBraunLab/gopeaks to your working directory</li>  
  <li>Download the appropriate species/assembly genome index files to your working directory (see below for sources)</li>  
</ol>    
  
#### These steps will be performed any time new data comes in:  
<ol>  
  <li>Upload (or download) your .fastq.gz read files to their respective directories</li>  
  <li>Modify the .sh scripts as needed to provide correct file paths and run options</li>  
  <li>Run the .sh script  
    <ul>  
      <li>sbatch ngs_pipeline_RNASeq.sh</li>  
    </ul>  
    or  
    <ul>  
      <li>sbatch ngs_pipeline_CUTnTag.sh</li>  
    </ul>  
  </li>  
</ol>  
___  

## Software Required  
<a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a> <a href="https://multiqc.info/">MultiQC</a>  
<a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> <a href="https://daehwankimlab.github.io/hisat2/">Hisat2</a>  
<a href="https://broadinstitute.github.io/picard/">Picard</a> <a href="http://www.htslib.org/">Samtools</a>  
<a href="https://deeptools.readthedocs.io/en/develop/">Deeptools</a>  
<a href="https://github.com/FredHutch/SEACR">SEACR</a> <a href="https://github.com/macs3-project/MACS">MACS</a> <a href="https://github.com/maxsonBraunLab/gopeaks">GoPeaks</a>  
  
On ComputeCanada (or other HPC platform), many of these programs are already available and you just have to load them. You can use the cc_modules_to_load.txt file shown below as a reference.  
  
### Loading software modules on ComputeCanada  
For python environment:  
>module load python>=3.8.2 scipy-stack/2021a  
  
For NGS tools:  
>module load samtools>=1.11 r>=4.2.1 bowtie2>=2.4.1 fastqc>=0.11.9 hisat2  

For featureCounts:  
>module load subread  

Meme Suite for motifs:  
>module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3  
>module load meme/5.4.1  

Extra tools converting filetypes e.g. bigwigtobedgraph:  
>module load nixpkgs/16.09  gcc/7.3.0  
>module load kentutils/20180716  

Bedtools:  
>module load bedtools/2.30.0  

Additionally, create a virtualenv using the provided cc_requirements.txt file if used on ComputeCanada, this environment must be loaded to run the script.  
>virutalenv ngsENV  
>source ngsENV/bin/activate  
>pip install -r cc_requirements.txt  
  
Example .sh scripts are also provided as templates for running on ComputeCanada or HPC platforms.  

___  
## Reference Genome Files  
Bowtie2 <a href="https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml">genome index files</a>  for alignment (on right under Indexes). Alternatively, you can download from <a href="https://support.illumina.com/sequencing/sequencing_software/igenome.html">iGenomes.</a>  
Hisat2 <a href="https://daehwankimlab.github.io/hisat2/download/">genome index files</a> for RNA-Seq alignment.  

### Blacklist Regions  
If blacklisted regions wish to be removed in bamCoverage, you can find these files <a href="https://github.com/Boyle-Lab/Blacklist">here.</a>  
___  
## Usage  
File name convention to run properly should follow hyphens "-" only between words and an underscore "_" before R1/R2 in filenames.  

&emsp;e.g.:  
&emsp;&emsp;<strong>CORRECT:</strong> one-file-of-reads_R1.fastq.gz  
&emsp;&emsp;<strong>NOT CORRECT:</strong> one_file-of_reads-R1.fastq.gz  
___  
#### 1. Prepare files  
Place .fastq.gz files for a sample replicate into a folder.  

&emsp;e.g. paired-end reads:  
&emsp;(Note: one paired-end replicate will have an "_R1" and "_R2" file)  
&emsp;&emsp;RAW_READS/  
&emsp;&emsp;&emsp;read_R1.fastq.gz  
&emsp;&emsp;&emsp;read_R2.fastq.gz  

&emsp;e.g. single-end reads:  
&emsp;(Note: only put multiple single-end replicates together if --merge will be used)  
&emsp;&emsp;RAW_READS/  
&emsp;&emsp;&emsp;read_1.fastq.gz  
&emsp;&emsp;&emsp;read_2.fastq.gz  
&emsp;&emsp;&emsp;read_3.fastq.gz  
&emsp;&emsp;&emsp;read_4.fastq.gz 
  
#### 2. Run command-line execution (see default options)    
&emsp;e.g. RNA-Seq:  
> python ngs_processing_pipeline.py -logfile RAW_READS.log -reads RAW_READS/ --technique rnaseq --species Mus --length 100 --reads_type paired --genome_index hisat_genomes_index/UCSC/mm10/genome --no_spikein -qctrim -outdir RAW_READS/  

&emsp;e.g. Other:  
> python ngs_processing_pipeline.py -logfile RAW_READS.log -reads RAW_READS/ --species Mus --length 100 --reads_type paired --genome_index mm10/Bowtie2Index/genome --no_spikein -qctrim -outdir RAW_READS/   

&emsp;e.g. Merge single-end replicates with "--merge" option:  
> python ngs_processing_pipeline.py -logfile RAW_READS.log -reads RAW_READS/ --species Mus --length 100 --reads_type single --genome_index mm10/Bowtie2Index/genome --no_spikein -qctrim -outdir RAW_READS/ --merge  
___ 
### Important output files  
&emsp;logs/  
&emsp;&emsp;&emsp;  ...log (monitor the progress of the script and troubleshoot problems)  
  
&emsp;Analysis_Results/QC_Rawreads/  
&emsp;&emsp;&emsp;  ...html (quality check raw reads and modify input options/re-run if required)  
  
&emsp;Analysis_Results/Peaks/  
&emsp;&emsp;&emsp; ...stringent.bed  
&emsp;&emsp;&emsp; ...peaks.narrowPeak (peaks files identifying enriched regions, useful in downstream analysis)  
&emsp;&emsp;&emsp; ...gopeaks_peaks.bed   
  
&emsp;Analysis_Results/Peaks/  
&emsp;&emsp;&emsp; ..._summits.bed (peak summits from MACS, useful in downstream analysis)  
  
&emsp;Analysis_Results/Normalized_and_Unnormalized_BigWigs/Normalized/  
&emsp;&emsp;&emsp; ...bw (normalized bigwigs for viewing coverage in genome browsers)  
  
&emsp;All_output/Processed_reads/  
&emsp;&emsp;&emsp; ...bam  
&emsp;&emsp;&emsp; ...bai  (alignment+index files (should always be together), required for many analysis tools)  
___
### Additional Resources  
<a href="https://www.nature.com/articles/s41467-019-09982-5">CUT&Tag</a>  
<a href="https://yezhengstat.github.io/CUTTag_tutorial/">Pipeline example</a>  
<a href="https://learn.gencore.bio.nyu.edu/">NGS Analysis</a>  
<a href="https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html">bamCoverage</a>  
<a href="https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_fastqc_assessment.html">FastQC Guide</a>  
<a href="https://www.genome.ucsc.edu/FAQ/FAQformat.html">File formats</a> and <a href="https://www.genome.ucsc.edu/FAQ/">UCSC Browser Help</a>  
<a href="https://igv.org/">Integrative Genomics Browser</a>  
<a href="https://docs.alliancecan.ca/wiki/Technical_documentation">Digital Research Alliance of Canada (Compute Canada) Technical Documentation</a>  
<a href="https://hgdownload.soe.ucsc.edu/admin/exe/">Binary downloads</a> for <a href="https://github.com/ucscGenomeBrowser/kent">kentutils</a> for converting file formats  
<a href="https://github.com/IARCbioinfo/BAM-tricks">BAM tricks</a>  

#### Genome Sequences  
<a href="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/">UCSC hg38 genome sequence files</a>  
<a href="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/">UCSC hg19 genome sequence files</a>  
<a href="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/">UCSC mm10 genome sequence files</a>  
___
### Pipeline  
<img src="ngs_pipeline.png" alt="Pipeline Process">    

___     
