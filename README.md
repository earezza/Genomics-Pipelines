# Genomics Pipelines  
### Description:  
Automated pipelines, scripts, and command-line operations to perform data processing, analysis, and visualizations of genomics data. 

- Process NGS .fastq.gz raw reads using an automated script to produce .bam, .bw, and .bed files. Executed on a high-performance computing platform.  
- Perform downstream analysis such as differential gene expression, gene annotations, and motif discovery.  
- Analyze single-cell RNA-Seq data by clustering, differential gene expression, and trajectory inference.  
  
### Prerequisites:    
- Linux terminal/bash/shell and <a href="https://ubuntu.com/tutorials/command-line-for-beginners#1-overview">command-line interface usage</a>  
- <a href="https://www.python.org/">Python</a> and <a href="https://docs.conda.io/en/latest/">conda</a>/<a href="https://virtualenv.pypa.io/en/latest/">virtualenv</a> package environment management  
- <a href="https://www.r-project.org/">R</a> and <a href="https://posit.co/downloads/">RStudio</a>  
- <a href="https://docs.jupyter.org/en/latest/">Jupyter Notebook</a>  
- <a href="https://slurm.schedmd.com/overview.html">Slurm</a> job scheduling  
  
### Resources  
NGS analysis NYU resources https://learn.gencore.bio.nyu.edu/  
Galaxy tools, training, tutorials https://usegalaxy.org/  
Bioinformatics Training HBC https://hbctraining.github.io/main/  

___  
 
#### A typical workflow for these scripts/tools may be:  
- For single-cell RNA-Seq data  
  1. Analyze-SingleCell-RNA-Seq/<a href="Analyze-SingleCell-RNA-Seq/scRNAseq_Notebook.ipynb">**scRNAseq_Notebook.ipynb**</a> or <a href="Analyze-SingleCell-RNA-Seq/scRNAseq_main.py">**scRNAseq_main.py**</a> (run steps manually as needed while exploring data)  
- For bulk RNA-Seq data  

  1. Process-NGS-RawReads/<a href="Process-NGS-RawReads/ngs_processing_pipeline.py">**ngs_processing_pipeline.py**</a> (with option -technique rnaseq)  
  2. DownstreamAnalysis/DifferentialGeneExpression/<a href="DownstreamAnalysis/DifferentialGeneExpression/RNA-Seq_Differential_Gene_Expression_Notebook.ipynb">**RNA-Seq_Differential_Gene_Expression_Notebook.ipynb**</a> or <a href="DownstreamAnalysis/DifferentialGeneExpression/analyze_rnaseq_degs.R">**analyze_rnaseq_degs.R**</a>  
  
- For ChIP-Seq, CUT&Tag, ATAC-Seq, etc... data  
  1. Process-NGS-RawReads/<a href="Process-NGS-RawReads/ngs_processing_pipeline.py>Process-NGS-RawReads/ngs_processing_pipeline.py">**ngs_processing_pipeline.py**</a> (with option -technique chipseq (or whichever data is from))  
  2. DownstreamAnalysis/DifferentialGeneExpression/<a href="DownstreamAnalysis/DifferentialGeneExpression/Peaks_Differential_Gene_Expression_Notebook.ipynb">**Peaks_Differential_Gene_Expression_Notebook.ipynb**</a> or <a href="DownstreamAnalysis/DifferentialGeneExpression/analyze_peaks_degs.R">**analyze_peaks_degs.R**</a>
  3. DownstreamAnalysis/MotifAnalysis/...
