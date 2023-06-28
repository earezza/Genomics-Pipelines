# This script requires a count matrix file produced by featureCounts and/or htseq-count
# and a sample info file to perform pairwise differential gene expression analysis using DESeq2.
# Additionally provides GO and KEGG annotations from DEG.

# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(clusterProfiler)
  library(ggplot2)
  library(EnhancedVolcano)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  #library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
  #library(ensembldb)
  #library(EnsDb.Hsapiens.v86)
  #library(EnsDb.Mmusculus.v79)
  #library(EnsDb.Rnorvegicus.v79)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(RColorBrewer)
  library(pheatmap)
  library(grid)
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
  make_option(c("-o", "--organism"), type="character", default="mouse", help="Organism to annotate genes (mouse (mm10) or human (hg38) or rat (rn6))", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="DEG_Analysis/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-f", "--filter"), action="store_true", type="logical", default=FALSE, help="Flag to filter read counts (removes genes < min_count from raw matrix, then removes genes < min_baseMean after normalizing", metavar="character"),
  make_option(c("--min_count"), type="integer", default=1, help="Minimum counts to include if filtering", metavar="integer"),
  make_option(c("--min_basemean"), type="double", default=10, help="Minimum baseMean of normalized counts to include if filtering", metavar="integer"),
  make_option(c("--lfc"), type="double", default=0.6, help="Magnitude of log2foldchange to define significant up/down regulation of genes", metavar="integer")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cat("Command-line options:\n")
for (i in which(names(opt) != "help")) {
  cat(names(opt)[i], '=', paste(opt)[i], "\n")
}
cat("log2FC of", opt$lfc, "equates to FC of", 2^opt$lfc)

# ========= SETUP VARIABLES =========
# Define references for GO/KEGG annotations
if (tolower(opt$organism) == "mouse"){
  annoDb <- "org.Mm.eg.db"
  keggOrg <- "mmu"
} else if (tolower(opt$organism) == "human"){
  annoDb <- "org.Hs.eg.db"
  keggOrg <- "hsa"
} else if (tolower(opt$organism) == "rat"){
  annoDb <- "org.Rn.eg.db"
  keggOrg <- "rno"
} else{
  stop("Invalid choice of organism")
}

# Define Functions
make_dotplot <- function(df, title="", ylabel="Description", colour="#56B1F7", n=15){
  df$ycolour <- "black"
  if ("ONTOLOGY" %in% colnames(df)){
    df$Description <- paste(df$ONTOLOGY, df$Description, sep=' - ')
    df$ycolour <- ifelse(grepl("BP -", df$Description), 'blue', df$ycolour)
    df$ycolour <- ifelse(grepl("CC -", df$Description), 'red', df$ycolour)
    df$ycolour <- ifelse(grepl("MF -", df$Description), 'darkgreen', df$ycolour)
  }
  df <- df[order(df$p.adjust, decreasing=FALSE),]
  plt <- ggplot() +
    geom_point(data=head(df, n=n),
               aes(x = -log(p.adjust), 
                   y = reorder(Description, -p.adjust), 
                   colour = Count, 
                   size = unname(unlist(sapply(GeneRatio, function(x) eval(parse(text=x)))))*100,
               ),
    ) + 
    theme_linedraw() +
    theme(axis.text.y = element_text(colour=rev(head(df$ycolour, n=n)))) +
    scale_color_gradient(low = "black", high = colour) +
    ggtitle(title)
  plt$labels$x <- "-log(p.adjust)"
  plt$labels$y <- ylabel
  plt$labels$size <- "GenePercentage"
  plt$labels$colour <- "GeneCount"
  return(plt)
}

make_volcanoplot <- function(res, condition1, condition2, lfc_threshold, padj_threshold){
  df <- as.data.frame(res)
  
  keyvals <- ifelse(df$log2FoldChange <= (0-lfc_threshold), 'darkred', 'black')
  keyvals <- ifelse(df$log2FoldChange <= (0-lfc_threshold) & df$padj <= padj_threshold, 'red', keyvals)
  keyvals <- ifelse(df$log2FoldChange >= lfc_threshold, 'darkgreen', keyvals)
  keyvals <- ifelse(df$log2FoldChange >= lfc_threshold & df$padj <= padj_threshold, 'green', keyvals)
  
  keyvals[is.na(keyvals)] <- 'black'
  names(keyvals)[keyvals == 'darkgreen'] <- condition1
  names(keyvals)[keyvals == 'green'] <- paste(condition1, '(significant)', sep=' ')
  names(keyvals)[keyvals == 'black'] <- 'Not DEG'
  names(keyvals)[keyvals == 'darkred'] <- condition2
  names(keyvals)[keyvals == 'red'] <- paste(condition2, '(significant)', sep=' ')
  
  plt <- EnhancedVolcano(df, 
                  x = 'log2FoldChange', 
                  y = 'padj', 
                  xlab = "Log2FoldChange",
                  ylab = "-Log(adjusted p-value)",
                  caption = paste0("Total = ", nrow(df), " genes\n", 
                                   condition1, " = ", nrow(subset(df, log2FoldChange >= lfc_threshold)), ' DEGs\n',
                                   condition2, " = ", nrow(subset(df, log2FoldChange <= (0-lfc_threshold))), ' DEGS'
                                   ),
                  colCustom = keyvals,
                  lab = rownames(df),
                  title = paste(condition1, ' vs ', condition2, sep=''),
                  subtitle = 'DESeq2 Results',
                  pCutoff = padj_threshold,
                  FCcutoff = lfc_threshold,
                  cutoffLineType = 'dashed',
                  legendLabels=c('Not sig.','Log2FC','p-value','p-value & Log2FC'),
                  labSize = 4,
                  legendPosition = 'top',
                  legendLabSize = 8,
                  )
  return(plt)
}

make_heatmapplot <- function(res, condition1, condition2, n = 50){
  df <- as.data.frame(head(res, n=n))
  df$Condition <- ifelse(df$log2FoldChange > 0, condition1, condition2)
  plt <- ggplot(df, aes(x = Condition, 
                 y = rownames(df), 
                 fill = log2FoldChange),
         ) +
    geom_tile() +
    ggtitle(paste("Top", n, "DEGs", sep=' '), subtitle='') +
    xlab("Condition") +
    ylab("Gene") 
  return(plt)
  
  
}

make_pheatmapplot <- function(anno, res, anno_type="GO", organism='mouse', heat_colour="PiYG", num_terms=25, num_genes=50, lfc=0.6, dendro=TRUE, sort_genes=TRUE, title="", xlabel="Gene", ylabel="Term"){
  
  # colour should be "Reds", "Greens", "Blues", or "PiYG"
  if ("ONTOLOGY" %in% colnames(anno)){
    anno$Description <- paste(anno$ONTOLOGY, anno$Description, sep=' - ')
  }
  
  # Take top n terms (most significant, already sorted by padj)
  df <- head(anno[order(anno$p.adjust, decreasing=FALSE), ], n=num_terms)
  
  # Create dataframe (matrix) of annotation terms vs genes with gene's associated log2FoldChange
  d <- data.frame()
  for (a in df$Description){
    gene_group <- strsplit(df[df$Description == a, ]$geneID, '/')[[1]]
    # For KEGG to convert EntrezID to gene Symbol
    if (anno_type == "KEGG"){
      if (tolower(organism) == "human"){
        gene_group <- mapIds(org.Hs.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }else if (tolower(organism) == "mouse"){
        gene_group <- mapIds(org.Mm.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }else if (tolower(organism) == "rat"){
        gene_group <- mapIds(org.Rn.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }
    }
    d[gene_group, a] <- res[gene_group, ]$log2FoldChange
  }
  # Sort by genes instead of by term (i.e. number of times gene found in all top terms)
  if (sort_genes){
    d <- d[names(sort(rowSums(is.na(d)))), ]
    xlabel <- paste(xlabel, "(Most Frequent in Top Terms)")
  }
  # Set NA to 0
  d[is.na(d)] <- 0
  
  range_max <- round(max(apply(d, 2, max)))
  range_min <- round(min(apply(d, 2, min)))
  #breaks <- seq( -round(2^lfc), round(2^lfc), length.out = 101)
  breaks <- seq( -abs((max(range_min, -round(lfc)))), min(range_max, round(lfc)), length.out = 101)
  #breaks <- seq( -abs((max(range_min, -round(2^lfc)))), min(range_max, round(2^lfc)), length.out = 101)
  if (range_max == 0){
    color <- rev(colorRampPalette(brewer.pal(n = 11, name = heat_colour))(101))
  }else{
    color <- colorRampPalette(brewer.pal(n = 11, name = heat_colour))(101)
  }
  
  # Plot
  setHook("grid.newpage", function() pushViewport(viewport(x=0,y=0.05,width=0.95, height=0.95, name="vp", just=c("left","bottom"))), action="prepend")
  pheatmap(t(head(d, n=num_genes)), 
           border_color = "grey90",
           color = color, # "Reds, Greens, Blues, RdYlGn for DEGs
           fontsize_row = 5,
           fontsize_col = 5,
           na_col = "white",
           breaks = breaks,
           cluster_rows = dendro,
           cluster_cols = dendro,
           main = title,
  ) 
  setHook("grid.newpage", NULL, "replace")
  grid.text("log2FoldChange", x=0.95, y=0.875, gp=gpar(fontsize=8))
  grid.text(xlabel, y=0, gp=gpar(fontsize=14))
  grid.text(ylabel, x=1, y=0.35,  rot=270, gp=gpar(fontsize=14))
  plt <- grid.grab()
  
  return(plt)
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
  cat("\n==========", comparison, "==========\n")
  
  # Isolate matrix for relevant samples
  samples1 <- rownames(subset(sampleinfo, Condition == combs[[c]][1]))
  samples2 <- rownames(subset(sampleinfo, Condition == combs[[c]][2]))
  mtx <- count_mtx[, c(samples1, samples2)]
  mtx_info <- subset(sampleinfo, rownames(sampleinfo) == samples1 | rownames(sampleinfo) == samples2)
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
    cat(quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), "\n")
    
    # Remove 0 counts (sum of all samples' counts) to avoid 0-bias and disregard irrelevant genes
    keep <- rowSums(counts(dds)) >= opt$min_count
    dds <- dds[keep,]
    cat("Percentiles after removing genes with sum(counts) <", opt$min_count, "\n")
    quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    
    # # Get relevant percentiles
    # quant <- quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    # # Remove genes with <= 0.1 percentile counts
    # cat("Removing genes with sum(counts) <=", quant[[2]])
    # keep <- rowSums(counts(dds)) > quant[[2]]
    # dds <- dds[keep,]
    
    # Keep genes only where there are at least 3 samples having at least that count (removes outliers)
    # cat("Removing genes with sum(counts) <=", quant[[2]], "in at least", round(length(colnames(dds))/4), "samples\n")
    # keep <- rowSums(counts(dds) > quant[[2]]) >= round(length(colnames(dds))/4)
    # dds <- dds[keep,]
    
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
  
  # Filter out normalized low-expressed genes (bottom 10th percentile) by baseMean  REMOVE ALL < 10 baseMean
  if (opt$filter){
    # cat("Genes and samples before filtering normalized genes:\n", dim(res), "\n")
    # quant <- quantile(res$baseMean, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    # cat("Removing genes with baseMean expression in bottom 10th percentile <=", quant[[2]], "\n")
    # res <- res[res$baseMean >= quant[[2]], ]
    # cat("Genes and samples after filtering normalized genes:\n", dim(res), "\n")
    
    # Remove any genes below minimum baseMean
    cat("Genes and samples before filtering normalized genes:\n", dim(res), "\n")
    cat("Removing genes with baseMean expression <", opt$min_basemean, "\n")
    res <- res[res$baseMean >= opt$min_basemean, ]
    cat("Genes and samples after filtering normalized genes:\n", dim(res), "\n")
  }
  
  # Generate plots
  png(paste(output_prefix, 'MAplot.png', sep=''))
  plotMA(res, ylim=c(-2,2), ylab='log2foldchange', alpha=0.1)
  title(main = comparison, sub="blue if adjusted pvalue < 0.1", cex.sub=0.8)
  abline(h=c((0-opt$lfc), opt$lfc), col="red", lwd=2)
  invisible(capture.output(dev.off()))
  
  resLFC <- lfcShrink(dds, coef=resultsNames(dds)[-1], type="apeglm")
  png(paste(output_prefix, 'shrunkMAplot.png', sep=''))
  plotMA(resLFC, ylim=c(-2,2), ylab='log2foldchange', alpha=0.1)
  title(main = comparison, sub="blue if adjusted pvalue < 0.1", cex.sub=0.8)
  abline(h=c((0-opt$lfc), opt$lfc), col="red", lwd=2)
  invisible(capture.output(dev.off()))
  
  plt <- make_volcanoplot(res, combs[[c]][1], combs[[c]][2], opt$lfc, 0.1)
  invisible(capture.output(ggsave(filename=paste(output_prefix, 'Volcanoplot.png', sep=''), plot=plt, dpi=320)))
  
  #png(paste(output_prefix, 'res_genecount_minpadj.png', sep=''))
  #plotCounts(dds, gene=which.min(res$padj), intgroup="Condition")
  #invisible(capture.output(dev.off()))
  
  #hist(res$padj[res$log2FoldChange >= 0.6], breaks = 0:20/20, col = "darkgreen", border = "white")
  #hist(res$padj[res$log2FoldChange <= -0.6], breaks = 0:20/20, col = "darkred", border = "white")
  
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
  
  # Heatmap
  plt <- make_heatmapplot(res_changed_sorted, combs[[c]][1], combs[[c]][2], n=40)
  invisible(capture.output(ggsave(filename=paste(output_prefix, 'Heatmapplot.png', sep=''), plot=plt, dpi=320)))
  
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
    }else if (tolower(opt$organism) == "rat"){
      genes_entrez[[n]] <- mapIds(org.Rn.eg.db, keys = genes[[n]], column = "ENTREZID", keytype = "SYMBOL")
    }
  }
  
  # ============= GO and KEGG Annotations ==============
  for (n in names(genes)){
    
    if (n == combs[[c]][1]){
      colour <- 'green'
      heat_colour <- "Greens"
    }else if (n== combs[[c]][2]){
      colour <- 'red'
      heat_colour <- "Reds"
    }else{
      colour <- "#56B1F7"
      heat_colour <- "PiYG"
    }
    
    # GO Annotation
    tryCatch(
      {
        for (ont in c('ALL', 'CC', 'MF', 'BP')){
          
          compGO <- compareCluster(geneCluster=genes[n],
                                   OrgDb=annoDb,
                                   fun="enrichGO",
                                   keyType="SYMBOL", #ENSEMBL
                                   pvalueCutoff=0.05,
                                   pAdjustMethod="BH",
                                   readable=TRUE,
                                   ont=ont
          ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          # Only output if there's results
          if ((!is.null(compGO)) & (dim(compGO@compareClusterResult)[1] > 0)){
            compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
            #plt <- dotplot(compGO, showCategory = 8, title = paste("GO -", n, sep=""))
            plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", n, sep=""), ylabel="GO Term", colour=colour, n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, 'GO_', ont, '_', n, '.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            # Write annotations to csv
            write.table(as.data.frame(compGO), file=paste(output_prefix, 'GO_', ont, '_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
            
            plt <- make_pheatmapplot(compGO@compareClusterResult, res, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE, title=paste("GO (", ont, ") - ", n, sep=""), ylabel="GO Term")
            invisible(capture.output(ggsave(filename=paste(output_prefix, 'GO_', ont, '_', n, '_pheatmap_by_gene.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            plt <- make_pheatmapplot(compGO@compareClusterResult, res, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE, title=paste("GO (", ont, ") - ", n, sep=""), ylabel="GO Term")
            invisible(capture.output(ggsave(filename=paste(output_prefix, 'GO_', ont, '_', n, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
          }
          
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
        if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
          #plt <- dotplot(compKEGG, showCategory = 8, title = paste("KEGG -", n, sep=""))
          plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', n, sep=""), ylabel="KEGG Category", colour=colour, n=15)
          invisible(capture.output(ggsave(filename=paste(output_prefix, 'KEGG_annotation_', n, '.png', sep=''), plot=plt, dpi=320)))
          remove(plt)
          # Write annotations to csv
          write.table(as.data.frame(compKEGG), file=paste(output_prefix, 'KEGG_annotation_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
          
          plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", organism=opt$organism, title=paste('KEGG - ', n, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE, ylabel="KEGG Category")
          invisible(capture.output(ggsave(filename=paste(output_prefix, 'KEGG_annotation_', n, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
          remove(plt)
          plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", organism=opt$organism, title=paste('KEGG - ', n, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE, ylabel="KEGG Category")
          invisible(capture.output(ggsave(filename=paste(output_prefix, 'KEGG_annotation_', n, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
          remove(plt)
        
        }
      },error = function(e)
      {
        message(e)
      }
    ) 
    
    # GSEA
    tryCatch(
      {
        res <- res[order(res$log2FoldChange, decreasing=TRUE),]
        original_gene_list <- res[genes[[n]],]$log2FoldChange
        names(original_gene_list) <- rownames(res[genes[[n]],])
        gene_list <- na.omit(original_gene_list)
        gene_list <- gene_list[order(gene_list, decreasing=TRUE)]
        
        for (ont in c('ALL', 'CC', 'MF', 'BP')){
          
          gsea <- gseGO(geneList=gene_list, 
                        ont = ont, 
                        keyType = "SYMBOL", 
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = annoDb, 
                        pAdjustMethod = "BH"
          )
          if ((!is.null(gsea)) & (dim(gsea@result)[1] > 0)){
            gsea@result <- gsea@result[order(gsea@result$p.adjust, decreasing=FALSE),] # Sort by most signiicant
            gsea@result <- head(gsea@result, n=20) # retain only top to plot
            plt <- dotplot(gsea, showCategory = 20, title = paste("GSEA (", ont, ") - ", n, sep=""))
            df <- plt$data
            plt <- make_dotplot(df, title=paste("GSEA (", ont, ") - ", n, sep=""), ylabel="GSEA", colour=colour, n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, 'GSEA_', ont, '_', n, '.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            
            # Write annotations to csv
            write.table(as.data.frame(gsea$result), file=paste(output_prefix, 'GSEA_', ont, '_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
            
            plt <- make_pheatmapplot(df, res, title=paste("GSEA (", ont, ") - ", n, sep=""), ylabel="GSEA", heat_colour=heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE)
            invisible(capture.output(ggsave(filename=paste(output_prefix, 'GSEA_', ont, '_', n, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            plt <- make_pheatmapplot(df, res, title=paste("GSEA (", ont, ") - ", n, sep=""), ylabel="GSEA", heat_colour=heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE)
            invisible(capture.output(ggsave(filename=paste(output_prefix, 'GSEA_', ont, '_', n, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
          }
          
        }
      },error = function(e)
      {
        message(e)
      }
    ) 
  }
  
}
