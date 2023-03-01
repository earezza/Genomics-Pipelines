
## Downstream Analysis: Differential Gene Expression  

These pipelines and tools are performed in R. Packages which may be required can be found in the providedd Renv_packages.txt file.  

### For RNA-Seq  
<a href="https://bioconductor.org/packages/release/bioc/html/Rsubread.html">Rsubread</a> to produce a count matrix from provided .bam files using featureCounts (if matrix not already created).   
<a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">DESeq2</a> to process the count matrix and produce differential gene expression statistics.   
<a href="https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html">clusterProfiler</a> to annotate genes (GO, GSEA, KEGG) using database packages for mappings.   
  
### For CUT&Tag, ChIP-Seq, etc...  
In addition to the above packages,  
<a href="https://bioconductor.org/packages/release/bioc/html/DiffBind.html">DiffBind</a> to compare peaksets, determine overlaps, and establish consensus peaksets.  

___  
### Additional Resources  
<a href="https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html">Normalization and DGE</a>  
<a href="https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html">RNA-Seq DGE example</a>  
<a href="https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html">RNA-Seq DGE example 2 with DESeq2</a>  
<a href="http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html">ChIPseeker examples</a>  
<a href="https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionV/lessons/12_annotation_functional_analysis.html">ChIPseeker examples 2</a>  
<a href="https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/09_data_visualization.html">deepTools visualizations examples</a>  
