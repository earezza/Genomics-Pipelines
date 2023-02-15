
## Downstream Analysis: Differential Gene Expression
Rsubread to produce a count matrix from provided .bam files using featureCounts (if matrix not already created). https://bioconductor.org/packages/release/bioc/html/Rsubread.html  
DESeq2 to process the count matrix and produce differential gene expression statistics. https://bioconductor.org/packages/release/bioc/html/DESeq2.html  
clusterProfilerto annotate genes (GO, GSEA, KEGG) using database packages for mappings. https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html  
