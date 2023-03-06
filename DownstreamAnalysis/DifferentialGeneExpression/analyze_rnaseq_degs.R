# This script requires a count matrix file produced by featureCounts and/or htseq-count
# and a sample info file to perform pairwise differential gene expression analysis using DESeq2.
# Additionally provides GO and KEGG annotations from DEG.

# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(clusterProfiler)
  library(ggplot2)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  #library(ensembldb)
  #library(EnsDb.Hsapiens.v86)
  #library(EnsDb.Mmusculus.v79)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(optparse)
}))

# If need to generate count matrix...Put all .bam and .bai files into a directory, e.g. "bams/"
#bam_files <- list.files(path="bams/", full.names=TRUE)[c(TRUE, FALSE)]
#bamcounts <- featureCounts(bam_files, annot.inbuilt="mm10", countMultiMappingReads=FALSE, ignoreDup=FALSE, isPairedEnd=TRUE, nthreads=4, verbose=TRUE)
#rownames(bamcounts$counts) <- mapIds(org.Mm.eg.db, keys = rownames(bamcounts$counts), column = "SYMBOL", keytype = "ENTREZID")
#bamcounts$counts <- bamcounts$counts[!(is.na(rownames(bamcounts$counts))), ]
#for (n in names(bamcounts)){
#  write.table(bamcounts[[n]], file=paste(getwd(), "/", n, ".csv", sep=""), sep=",", quote=F, col.names=NA)
#}

# ======= Get command-line optional arguments =======
option_list = list(
  make_option(c("-c", "--countsfile"), type="character", help="Count matrix file", metavar="character"),
  make_option(c("-s", "--sampleinfo"), type="character", help="File containing samples and their experimental conditions info (row names must match matrix file column names)", metavar="character"),
  make_option(c("-o", "--organism"), type="character", default="mouse", help="Organism to annotate genes (mouse or human)", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="DEG_Analysis/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-f", "--filter"), action="store_true", type="logical", default=FALSE, help="Flag to filter read counts (removes 0 counts, then removes genes in bottom 10th percentile of counts for those in < 3 samples)", metavar="character"),
  make_option(c("--min_count"), type="integer", default=1, help="Minimum counts to include if filtering", metavar="integer"),
  make_option(c("--lfc"), type="double", default=1.5, help="Magnitude of log2foldchange to define significant up/down regulation of genes", metavar="integer")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cat("Command-line options:\n")
for (i in which(names(opt) != "help")) {
  cat(names(opt)[i], '=', paste(opt)[i], "\n")
}

# ========= SETUP VARIABLES =========
# Define references for GO/KEGG annotations
if (tolower(opt$organism) == "mouse"){
  annoDb <- "org.Mm.eg.db"
  keggOrg <- "mmu"
} else if (tolower(opt$organism) == "human"){
  annoDb <- "org.Hs.eg.db"
  keggOrg <- "hsa"
} else{
  stop("Invalid choice of organism")
}

# =========== Load Input Files ============
count_mtx <- as.matrix(read.csv(opt$countsfile, sep=",", row.names=1, check.names=FALSE))
sampleinfo <- read.csv(opt$sampleinfo, row.names=1)

# Create output file directory and set as working directory
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}
setwd(opt$result_dir)
cat("Output files will be in", paste(getwd(), "/", sep=""), "\n")

# Get all pairwise combinations of conditions for comparing
combs <- as.data.frame(combn(unique(sampleinfo$Condition), 2))

# Iterate through condition combinations to compare DEG
for (c in colnames(combs)){
  
  comparison <- paste(combs[[c]][1], combs[[c]][2], sep="_vs_")
  
  # Isolate matrix for relevant samples
  samples1 <- rownames(sampleinfo[sampleinfo$Condition == combs[[c]][1], ])
  samples2 <- rownames(sampleinfo[sampleinfo$Condition == combs[[c]][2], ])
  mtx <- count_mtx[, c(samples1, samples2)]
  mtx_info <- sampleinfo[c(samples1, samples2), ]
  mtx_info$Condition <- factor(mtx_info$Condition)
  
  # Check that matrix columns are same order as sample info
  if (!all(rownames(mtx_info) == colnames(mtx))) {
    break
  }
  
  # Create DESeqDataSet object from input
  dds <- DESeqDataSetFromMatrix(countData = mtx,
                                colData = mtx_info,
                                design = ~ Condition)
  
  # Create directory for output files of current comparison
  if (!file.exists(comparison)) {
    dir.create(comparison)
  }
  
  # Set output filenames prefix
  output_prefix <- paste(getwd(), "/", comparison, "/", sep="")
  
  cat("Genes and samples:\n", dim(dds), "\n")
  
  # =========== Filter Count Matrix ============
  if (opt$filter) {
    
    # Get percentile ranges of counts
    cat("Percentiles of genes count sums over all samples:\n")
    quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    
    # Remove 0 counts to avoid 0-bias and disregard irrelevant genes
    keep <- rowSums(counts(dds)) >= opt$min_count
    dds <- dds[keep,]
    cat("Percentiles after removing genes with sum(counts) <", opt$min_count, "\n")
    quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    # Get relevant percentiles
    quant <- quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    # Remove genes with <= 0.1 percentile counts
    cat("Removing genes with sum(counts) <=", quant[[2]])
    keep <- rowSums(counts(dds)) > quant[[2]]
    dds <- dds[keep,]
    # Keep genes only where there are at least 3 samples having at least that count
    cat("Removing genes with sum(counts) <=", quant[[2]], "in at least", round(length(colnames(dds))/4), "samples\n")
    keep <- rowSums(counts(dds) > quant[[2]]) >= round(length(colnames(dds))/4)
    dds <- dds[keep,]
    
    cat("Genes and samples after filtering raw counts:\n", dim(dds), "\n")
  }
  
  # ========== Normalize and Compute DEG Stats =============
  # Create DESeq object
  dds <- DESeq(dds)
  
  # Write filtered count matrix to file
  write.table(counts(dds), file=paste(output_prefix, "count_mtx_filtered.csv", sep=""), sep=",", quote=F, col.names=NA)
  
  # Write normalized values to file
  normalized_count_mtx <- counts(dds, normalized=TRUE)
  write.table(normalized_count_mtx, file=paste(output_prefix, "count_mtx_normalized.csv", sep=""), sep=",", quote=F, col.names=NA)
  
  # Obtain DEG results and output to file
  res <- results(dds, contrast=c("Condition", combs[[c]][1], combs[[c]][2]))
  res[['FoldChange']] <- 2^abs(res[['log2FoldChange']])*(res[['log2FoldChange']]/abs(res[['log2FoldChange']]))
  write.table(res[order(res$log2FoldChange, decreasing=TRUE), ], file=paste(output_prefix, "DESeq2_FullResult_", comparison, ".csv", sep=""), sep=",", quote=F, col.names=NA)
  
  #resLFC <- lfcShrink(dds, coef=resultsNames(dds)[-1], type="apeglm")
  
  # Filter out normalized low-expressed genes (bottom 10th percentile) by baseMean
  cat("Genes and samples before filtering normalized genes:\n", dim(res), "\n")
  quant <- quantile(res$baseMean, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
  res <- res[res$baseMean >= quant[[2]], ]
  cat("Genes and samples after filtering normalized genes:\n", dim(res), "\n")
  
  # Generate plots
  png(paste(output_prefix, 'MAplot.png', sep=''))
  plotMA(res)
  abline(h=c((0-opt$lfc), opt$lfc), col="red", lwd=2)
  invisible(capture.output(dev.off()))
  
  #png(paste(output_prefix, 'resLFC_MAplot.png', sep=''))
  #plotMA(resLFC, ylim=c(-2,2))
  #invisible(capture.output(dev.off()))
  
  #png(paste(output_prefix, 'res_genecount_minpadj.png', sep=''))
  #plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")
  #invisible(capture.output(dev.off()))
  
  vsd <- vst(dds, blind=FALSE)
  png(paste(output_prefix, 'PCAplot.png', sep=''))
  plotPCA(vsd, intgroup=c("Condition"))
  invisible(capture.output(dev.off()))
  
  # ========== Organize DEGs ==============
  # Remove any NAs
  #res <- res[rowSums(is.na(res)) == 0, ]
  
  # Get significantly Upregulated (log2FoldChange > 0)
  res_up <- res[res$log2FoldChange >= opt$lfc ,]
  res_up_sorted <- res_up[order(res_up$log2FoldChange, decreasing=TRUE),]
  write.table(res_up_sorted, file=paste(output_prefix, "DESeq2_Result_", combs[[c]][1], ".csv", sep=""), sep=",", quote=F, col.names=NA)
  
  # Get significantly Downregulated (log2FoldChange < 0)
  res_down <- res[res$log2FoldChange <= (0 - opt$lfc) ,]
  res_down_sorted <- res_down[order(res_down$log2FoldChange, decreasing=TRUE),]
  write.table(res_down_sorted, file=paste(output_prefix, "DESeq2_Result_", combs[[c]][2], ".csv", sep=""), sep=",", quote=F, col.names=NA)
  
  # Get all significantly changed genes
  res_changed <- res[(res$log2FoldChange >= opt$lfc | res$log2FoldChange <= (0 - opt$lfc)),]
  res_changed_sorted <- res_changed[order(abs(res_changed$log2FoldChange), decreasing=TRUE),]
  write.table(res_changed_sorted, file=paste(output_prefix, "DESeq2_Result_DEGs", ".csv", sep=""), sep=",", quote=F, col.names=NA)
  
  genes <- list()
  genes[[combs[[c]][1]]] <- rownames(res_up_sorted)
  genes[[combs[[c]][2]]] <- rownames(res_down_sorted)
  genes[["DEGs"]] <- rownames(res_changed_sorted)
  
  genes_entrez <- list()
  for (n in names(genes)){
    if (tolower(opt$organism) == "human"){
      genes_entrez[[n]] <- mapIds(org.Hs.eg.db, keys = genes[[n]], column = "ENTREZID", keytype = "SYMBOL")
    }else if (tolower(opt$organism) == "mouse"){
      genes_entrez[[n]] <- mapIds(org.Mm.eg.db, keys = genes[[n]], column = "ENTREZID", keytype = "SYMBOL")
    }
  }
  
  # ============= GO and KEGG Annotations ==============
  for (n in names(genes)){
    
    # GO Annotation
    tryCatch(
      {
        compGO <- compareCluster(geneCluster=genes[n],
                                 OrgDb=annoDb,
                                 fun="enrichGO",
                                 keyType="SYMBOL", #ENSEMBL
                                 pvalueCutoff=0.05,
                                 pAdjustMethod="BH",
                                 readable=TRUE
        ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
        if (!is.null(compGO)){
          plt <- dotplot(compGO, showCategory = 5, title = paste("GO -", n, sep=""))
          invisible(capture.output(ggsave(filename=paste(output_prefix, 'GO_annotation_', n, '.png', sep=''), plot=plt, dpi=320)))
          # Write annotations to csv
          write.table(as.data.frame(compGO), file=paste(output_prefix, 'GO_annotation_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
        }
      },error = function(e)
      {
        message(e)
      }
    )
    
    # KEGG Annotation
    tryCatch(
      {
        # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
        compKEGG <- compareCluster(geneCluster=genes_entrez[n],
                                   fun="enrichKEGG",
                                   pvalueCutoff=0.05,
                                   pAdjustMethod="BH",
                                   organism=keggOrg
        ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
        if (!is.null(compKEGG)){
          plt <- dotplot(compKEGG, showCategory = 10, title = paste("KEGG -", n, sep=""))
          invisible(capture.output(ggsave(filename=paste(output_prefix, 'KEGG_annotation_', n, '.png', sep=''), plot=plt, dpi=320)))
          
          # Write annotations to csv
          write.table(as.data.frame(compKEGG), file=paste(output_prefix, 'KEGG_annotation_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
        }
      },error = function(e)
      {
        message(e)
      }
    ) 
    
    # GSEA
    tryCatch(
      {
        original_gene_list <- res[genes[[n]],]$log2FoldChange
        names(original_gene_list) <- rownames(res[genes[[n]],])
        gene_list <- na.omit(original_gene_list)
        
        gsea <- gseGO(geneList=gene_list, 
                      ont ="ALL", 
                      keyType = "SYMBOL", 
                      nPerm = 10000, 
                      minGSSize = 3, 
                      maxGSSize = 800, 
                      pvalueCutoff = 0.05, 
                      verbose = TRUE, 
                      OrgDb = annoDb, 
                      pAdjustMethod = "BH"
        )
        if (!is.null(gsea)){
          plt <- dotplot(gsea, showCategory = 10, title = paste("GSEA -", n, sep=""))
          invisible(capture.output(ggsave(filename=paste(output_prefix, 'GSEA_annotation_', n, '.png', sep=''), plot=plt, dpi=320)))
          
          # Write annotations to csv
          write.table(as.data.frame(compKEGG), file=paste(output_prefix, 'GSEA_annotation_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
        }
      },error = function(e)
      {
        message(e)
      }
    ) 
  }
  
}
