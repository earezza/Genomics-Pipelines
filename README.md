# CnT-Pipeline

## Software Required  
FASTQC https://www.bioinformatics.babraham.ac.uk/projects/fastqc/  
MULTIQC https://multiqc.info/  
Bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml  
Picard https://broadinstitute.github.io/picard/  
Samtools  http://www.htslib.org/  
Deeptools  https://deeptools.readthedocs.io/en/develop/  
SEACR https://github.com/FredHutch/SEACR  
MACS2 https://github.com/macs3-project/MACS/wiki/Install-macs2  
GoPeaks https://github.com/maxsonBraunLab/gopeaks  
  
### Resources  
CUT&Tag https://www.nature.com/articles/s41467-019-09982-5  
Pipeline example https://yezhengstat.github.io/CUTTag_tutorial/  
Additional https://learn.gencore.bio.nyu.edu/  
___
Start  
	  Reads_fastq.gz  
	  Reads.fastq.gz.md5  

1. md5checksum (ensure files not corrupted in transfers, copies, etc...)  

2. QC -> fastqc (+picard)	(quality control of reads)  
	  Reads_fastqc.html  
	  Reads_fastqc.zip  

3. CompileQCResults -> multiqc (searches directory for results and generates report into .html)  
	  Rawreads_QC.html  

4. AdapterTrim -> cutadapt (finds and removes adapter sequences, primers, poly-A tails, etc...from reads)  
	  Reads_Trimmed.fastq  

5. CompileResultsPostTrimQC -> multiqc  
	  postTrimming_QC.html  

6. MapBowtie2 -> bowtie2, picard, samtools (bowtie2 aligns reads to reference sequences, samtools creates .bam file and .bam.bai index file, picard sorts .bam file by coordinate)  
	  Reads.bam  
	  Reads.coordsorted.bam  
	  Reads.coordsorted.bam.bai  

7. CollectAlignmentStats -> picard, samtools (marks duplicate reads in .bam file and generates new with duplicates tagged)  
	  Reads.dupMarked.bam  
	  Reads_picard.dupMark.txt  

8. CompileResultsMap -> multiqc  
	  Alignment_results.html  

9. FilteringBamsPicardSamtools -> samtools, picard (skip alignments with mapping quality scoring < 10 −10 log10 Pr{mapping position is wrong}, remove duplicates)  
	  Reads.MappedPaired.MAPQ10.bam  
	  Reads.MappedPaired.MAPQ10.bam.bai  
	  Reads.MappedPaired.MAPQ10.NoDups.bam  	
	  Reads_picard.rmDup.txt  
	  Reads.MappedPaired.MAPQ10.NoDups.bam.bai  

10. CompileResultsFiltering -> multiqc  
	  filteringbamsStats.html  

11. GetBigwigsBamCoverage -> deeptools (bamCoverage) (convert .bam to .bw (or .bed))  
	  Reads_RPGC.bw  
	  Reads_wo.norm.bw  
	  Reads_wo.norm_wDups.bw  

12. Map2Spikein_Bowtie2 -> bowtie2, picard, samtools (alignment with spike-in index)  
	  Reads.bam  
	  Reads.coordsorted.bam  
	  Reads.coordsorted.bam.bai  

13. CollectSpikeAlignmentStats -> picard, samtools (mark duplicates)  
	  Reads.dupMarked.bam  
	  Reads_picard.dupMark.txt  

14. CompileResultsSpike -> multiqc  
	  Spike_alignment.html  

15. CalcNormFactors	(get scaling to normalize reads to spike-in)  
	  Spike_align_stats.csv  

16. GetNormBwsBdgsBamCoverage -> deeptools (get normalized reads as .bw with/without duplicates, get bedgraph for SEACR)  
	  Reads_Norm.bw  
	  Reads_Norm.bedgraph  
	  Reads_Norm_wDups.bw  

17. Peak_Calling -> SEACR (peaks by calling enriched regions, can use control and stringent or relaxed thresholds), MACS2, GoPeaks  
	  Reads.stringent.bed  
	  Reads.relaxed.bed  
	  Reads_peaks.narrowPeak  
	  Reads_gopeaks_peaks.bed  
	  

18. Cleanup  
End  
