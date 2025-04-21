
## Downstream Analysis: Differential Binding Sites and Gene Expression  

These pipelines and tools are performed in R. Packages which may be required can be found in the provided <strong>Renv_packages.txt</strong> file.  
  
R version 4.3.1  
or download my_packages.txt and perform the following from command-line with your newer R version:  
Rscript -e 'packages <- readLines("my_packages.txt"); install.packages(packages)'  
  
### For RNA-Seq (analyze_rnaseq_degs_<DESeq2 or edgeR>.R)  
#### Create count matrix with get_rnaseq_counts.R  
<a href="https://bioconductor.org/packages/release/bioc/html/Rsubread.html">Rsubread</a> to produce a count matrix from provided .bam files using featureCounts (if matrix not already created).   
#### Run analysis  
<a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">DESeq2</a> to process the count matrix and produce differential gene expression statistics by median of ratios.    
<a href="https://bioconductor.org/packages/release/bioc/html/edgeR.html">edgeR</a> to process the count matrix and produce differential gene expression statistics by trimmed mean of M values.    
<a href="https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html">clusterProfiler</a> and <a href="https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html">ChIPseeker</a> to annotate genes (GO, GSEA, KEGG) using database packages for mappings.   
  
### For CUT&Tag, ChIP-Seq, etc... (analyze_peaks_degs.R)  
In addition to the above packages,  
<a href="https://bioconductor.org/packages/release/bioc/html/DiffBind.html">DiffBind</a> to compare peaksets, determine overlaps, and establish consensus peaksets.  

___  
### Additional Resources  
<a href="https://yulab-smu.top/biomedical-knowledge-mining-book/index.html">Biomedical Knowledge Mining using GOSemSim and clusterProfiler</a>  
<a href="https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html">DGE with DESeq2</a>  
<a href="https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html">Normalization and DGE</a>  
<a href="https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html">RNA-Seq DGE example</a>  
<a href="https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html">RNA-Seq DGE example 2 with DESeq2</a>  
<a href="http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html">ChIPseeker examples</a>  
<a href="https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html">ChIPseeker examples 2</a>  
<a href="https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/09_data_visualization.html">deepTools visualizations examples</a>  
