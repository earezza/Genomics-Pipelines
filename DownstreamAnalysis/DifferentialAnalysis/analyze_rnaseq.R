# This script requires a count matrix file produced by featureCounts and/or htseq-count
# and a sample info file to perform pairwise differential gene expression analysis using DESeq2.
# Additionally provides GO and KEGG annotations from DEG.

# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(DESeq2)
  library(edgeR)
  library(clusterProfiler)
  library(ggplot2)
  library(EnhancedVolcano)
  #library(ensembldb)
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
  make_option(c("-a", "--assembly"), type="character", default="mm10", help="Assembly to annotate genes/peaks (e.g. hg19, hg38, mm9, mm10, rn6)", metavar="character"),
  make_option(c("-d", "--database"), type="character", default="ucsc", help="Database reference for peaks gene annotations, ucsc (default) or ensembl", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="DEG_Analysis/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-m", "--method"), type="character", default="all", help="Method to analyze data, ('deseq2', 'edger', default: 'all')", metavar="character"), 
  make_option(c("-f", "--filter"), action="store_true", type="logical", default=FALSE, help="Flag to filter read counts (removes genes < min_count from raw matrix, then removes genes < min_baseMean after normalizing", metavar="character"),
  make_option(c("--min_count"), type="integer", default=1, help="Minimum counts to include if filtering", metavar="integer"),
  make_option(c("--min_basemean"), type="double", default=10, help="Minimum baseMean of normalized counts to include if filtering", metavar="integer"),
  make_option(c("--lfc"), type="double", default=0.585, help="Magnitude of log2foldchange to define significant up/down regulation of genes", metavar="integer"),
  make_option(c("--pvalue"), type="double", default=0.05, help="Significance threshold for DEGs (pvalue instead of p.adjusted for case where small sample set results in p.adjust=NA)", metavar="integer")
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
# Define Functions
load_annotation <- function(assembly, database){
  # ========= Get database references for annotations =========
  if (assembly == "mm10" | assembly == "mm9"){
    library(org.Mm.eg.db)
    annoDb <- "org.Mm.eg.db"
    keggOrg <- "mmu"
    if (database == "ucsc"){
      if (assembly == "mm10"){
        library(TxDb.Mmusculus.UCSC.mm10.knownGene)
        txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
      }else if (assembly == "mm9"){
        library(TxDb.Mmusculus.UCSC.mm9.knownGene)
        txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
      }
    }else if (database == "ensemble"){
      library(EnsDb.Mmusculus.v79)
      txdb <- EnsDb.Mmusculus.v79
      seqlevelsStyle(txdb) <- "UCSC" # format ensembl genes using UCSC style
    }
    else{
      stop("Invalid choice of annotation database")
    }
  } else if (assembly == "hg38" | assembly == "hg19"){
    library(org.Hs.eg.db)
    annoDb <- "org.Hs.eg.db"
    keggOrg <- "hsa"
    if (database == "ucsc"){
      if (assembly == "hg38"){
        library(TxDb.Hsapiens.UCSC.hg38.knownGene)
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
      }else if (assembly == "hg19"){
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
      }
    }else if (database == "ensemble"){
      library(EnsDb.Hsapiens.v86)
      txdb <- EnsDb.Hsapiens.v86
      seqlevelsStyle(txdb) <- "UCSC" # format ensembl genes using UCSC style
    } 
    else{
      stop("Invalid choice of annotation database")
    }
  } else if (assembly == "rn6"){
    library(TxDb.Rnorvegicus.UCSC.rn6.refGene)
    library(org.Rn.eg.db)
    annoDb <- "org.Rn.eg.db"
    keggOrg <- "rno"
    if (database == "ucsc"){
      txdb <- TxDb.Rnorvegicus.UCSC.rn6.refGene
    }else if (database == "ensemble"){
      library(EnsDb.Rnorvegicus.v79)
      txdb <- EnsDb.Rnorvegicus.v79
      seqlevelsStyle(txdb) <- "UCSC" # format ensembl genes using UCSC style
    } 
    else{
      stop("Invalid choice of annotation database")
    }
  } else{
    stop("Invalid choice of assembly")
  }
  
  anno_ref <- list("txdb"=txdb, "annoDb"=annoDb, "keggOrg"=keggOrg)
  
  return(anno_ref)
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
  
  # Re-format y-axis labels to not squish graph
  for (d in 1:length(df$Description)){
    i <- 1
    s <- df$Description[d]
    df$Description[d] <- ""
    while (i < length(strsplit(s, ' ')[[1]])){
      df$Description[d] <- paste(df$Description[d], paste(strsplit(s, ' ')[[1]][i:(i+5)], collapse = ' '), sep='\n')
      i <- (i+6)
    }
    df$Description[d] <- gsub(" NA", "", df$Description[d])
    df$Description[d] <- substring(df$Description[d], 2, nchar(df$Description[d]))
  }
  # Plot
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

make_volcanoplot <- function(res, condition1, condition2, lfc_threshold, padj_threshold, subtitle=''){
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
                         subtitle = subtitle,
                         pCutoff = padj_threshold,
                         FCcutoff = lfc_threshold,
                         cutoffLineType = 'dashed',
                         pointSize = 2.0,
                         legendLabels=c('Not sig.','Log2FC','p-value','p-value & Log2FC'),
                         labSize = 3,
                         legendPosition = 'top',
                         legendLabSize = 6,
  )
  return(plt)
}

make_heatmapplot <- function(res, condition1, condition2, n = 50, subtitle=''){
  df <- as.data.frame(head(res, n=n))
  df$Condition <- ifelse(df$log2FoldChange > 0, condition1, condition2)
  plt <- ggplot(df, aes(x = Condition, 
                        y = rownames(df), 
                        fill = log2FoldChange),
  ) +
    geom_tile() +
    ggtitle(paste("Top", n, "DEGs", sep=' '), subtitle=subtitle) +
    xlab("Condition") +
    ylab("Gene") 
  return(plt)
  
  
}

make_pheatmapplot <- function(anno, res, anno_type="GO", assembly='mm10', heat_colour="PiYG", num_terms=25, num_genes=50, lfc=0.6, dendro=TRUE, sort_genes=TRUE, title="", xlabel="Gene", ylabel="Term"){
  
  # colour should be "Reds", "Greens", "Blues", or "PiYG"
  if ("ONTOLOGY" %in% colnames(anno)){
    anno$Description <- paste(anno$ONTOLOGY, anno$Description, sep=' - ')
  }
  
  # Take top n terms (most significant, already sorted by padj)
  df <- head(anno[order(anno$p.adjust, decreasing=FALSE), ], n=num_terms)
  
  # Re-format y-axis labels to not squish graph
  for (d in 1:length(df$Description)){
    i <- 1
    s <- df$Description[d]
    df$Description[d] <- ""
    while (i < length(strsplit(s, ' ')[[1]])){
      df$Description[d] <- paste(df$Description[d], paste(strsplit(s, ' ')[[1]][i:(i+5)], collapse = ' '), sep='\n')
      i <- (i+6)
    }
    df$Description[d] <- gsub(" NA", "", df$Description[d])
    df$Description[d] <- substring(df$Description[d], 2, nchar(df$Description[d]))
  }
  
  # Create dataframe (matrix) of annotation terms vs genes with gene's associated log2FoldChange
  d <- data.frame()
  for (a in df$Description){
    gene_group <- strsplit(df[df$Description == a, ]$geneID, '/')[[1]]
    # For KEGG to convert EntrezID to gene Symbol
    if (anno_type == "KEGG"){
      if (assembly == "hg19" | assembly == "hg38"){
        gene_group <- mapIds(org.Hs.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }else if (assembly == "mm9" | assembly == "mm10"){
        gene_group <- mapIds(org.Mm.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }else if (assembly == "rn6"){
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

anno_ref <- load_annotation(opt$assembly, opt$database)

# =========== Load Input Files ============
count_mtx <- as.matrix(read.csv(opt$countsfile, sep=",", row.names=1, check.names=FALSE))
sampleinfo <- read.csv(opt$sampleinfo, row.names=1)

count_mtx <- count_mtx[, rownames(sampleinfo)]

# Create output file directory and set as working directory
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}
#setwd(opt$result_dir)
cat("Output files will be in", opt$result_dir, "\n")


change_dirs <- function(res_dir, method, subfolder){
  if ((method == 'DESeq2')){
    method_dir <- paste(res_dir, 'DESeq2/', sep='')
    if (!file.exists(paste(method_dir, subfolder, '/', sep=''))) {
      dir.create(paste(method_dir, subfolder, '/', sep=''))
    }
  }else if ((method == 'edgeR')){
    method_dir <- paste(res_dir, 'edgeR/', sep='')
    if (!file.exists(paste(method_dir, subfolder, '/', sep=''))) {
      dir.create(paste(method_dir, subfolder, '/', sep=''))
    }
  }else{
    method_dir <- paste(res_dir, method, '/', sep='')
    if (!file.exists(paste(method_dir, subfolder, '/', sep=''))) {
      dir.create(paste(method_dir, subfolder, '/', sep=''))
    }
  }
  return(paste(method_dir, subfolder, '/', sep=''))
}


# Get all pairwise combinations of conditions for comparing
combs <- as.data.frame(combn(unique(sampleinfo$Condition), 2))


# Perform global comparison (obtain PCA plot)
# DESeq2
if (opt$method == 'deseq2' | opt$method == 'all') {
  output_prefix <- change_dirs(opt$result_dir, 'DESeq2', '')
  dds <- DESeqDataSetFromMatrix(countData = count_mtx,
                                colData = sampleinfo,
                                design = ~ Condition)
  if (opt$filter) {
    # Get percentile ranges of counts
    cat("Filtering counts...\n")
    # Remove 0 counts (sum of all samples' counts) to avoid 0-bias and disregard irrelevant genes
    keep <- rowSums(counts(dds)) >= opt$min_count
    dds <- dds[keep,]
  }
  dds <- DESeq(dds)
  
  vsd <- vst(dds, blind=FALSE)
  plt <- plotPCA(vsd, intgroup=c("Condition")) + 
    geom_text(aes(label = rownames(vsd@colData)), nudge_y = 0.2, size = 3)
  invisible(capture.output(ggsave(paste(output_prefix, "PCAplot.png", sep=''), plot=plt)))
}

# edgeR
if (opt$method == 'edger' | opt$method == 'all') {
  output_prefix <- change_dirs(opt$result_dir, 'edgeR', '')
  edge <- DGEList(counts=count_mtx, group=sampleinfo$Condition, remove.zeros=FALSE)
  if (opt$filter) {
    cat("Filtering counts...\n")
    keep <- rowSums(edge$counts) >= opt$min_count
    edge <- edge[keep,]
  }
  edge <- calcNormFactors(edge, method="TMM")
  normalized_count_mtx_edge <- cpm(edge, normalized.lib.sizes=TRUE) 
  edge <- DGEList(counts = normalized_count_mtx_edge, group=sampleinfo$Condition, remove.zeros=FALSE, norm.factors = edge$samples$norm.factors)
  
  # Plot PCA
  png(paste(output_prefix, 'MDSplot.png', sep=''))
  plotMDS(edge, col=as.numeric(edge$samples$group), plot=TRUE)
  invisible(capture.output(dev.off()))
}

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
  if (opt$method == 'deseq2' | opt$method == 'all') {
    dds <- DESeqDataSetFromMatrix(countData = mtx,
                                  colData = mtx_info,
                                  design = ~ Condition)
    
    cat("DESeq2 genes and samples:\n", dim(dds), "\n")
  }
  
  # Create edgeR object
  if (opt$method == 'edger' | opt$method == 'all') {
    edge <- DGEList(counts=mtx, 
                    group=mtx_info$Condition, 
                    remove.zeros=FALSE)
    
    cat("EdgeR genes and samples:\n", dim(edge), "\n")
  }
  
  # =========== Filter Count Matrix ============
  if (opt$filter) {
    
    # Get percentile ranges of counts
    if (opt$method == 'deseq2' | opt$method == 'all') {
      cat("DESEq2 percentiles of genes count sums over all samples:\n")
      cat(quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), "\n")
      
      # Remove 0 counts (sum of all samples' counts) to avoid 0-bias and disregard irrelevant genes
      keep_dds <- rowSums(counts(dds)) >= opt$min_count
      dds <- dds[keep_dds,]
      cat("DESeq2 percentiles after removing genes with sum(counts) <", opt$min_count, "\n")
      cat(quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)))
      low_counts <- quantile(rowSums(counts(dds)), probs = c(0.10))
      
      # Remove bottom 10 percentile counts from remaining
      keep_dds <- rowSums(counts(dds)) > low_counts
      dds <- dds[keep_dds,]
      cat("DESeq2 percentiles after removing genes with sum(counts) <=", low_counts, "\n")
      cat(quantile(rowSums(counts(dds)), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)))
      
      cat("DESeq2 genes and samples after filtering raw counts:\n", dim(dds), "\n")
    }
    
    
    if (opt$method == 'edger' | opt$method == 'all') {
      cat("EdgeR percentiles of genes count sums over all samples:\n")
      cat(quantile(rowSums(edge$counts), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)), "\n")
    
      keep_edge <- rowSums(edge$counts) >= opt$min_count
      edge <- edge[keep_edge,]
      cat("EdgeR percentiles after removing genes with sum(counts) <", opt$min_count, "\n")
      cat(quantile(rowSums(edge$counts), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)))
      
      low_counts <- quantile(rowSums(edge$counts), probs = c(0.10))
      
      # Remove bottom 10 percentile counts from remaining
      keep_edge <- rowSums(edge$counts) > low_counts
      edge <- edge[keep_edge,]
      cat("DESeq2 percentiles after removing genes with sum(counts) <=", low_counts, "\n")
      cat(quantile(rowSums(edge$counts), probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)))
      cat("EdgeR genes and samples after filtering raw counts:\n", dim(edge), "\n")
    }
    
  }
  
  # ========== Normalize and Compute DEG Stats =============
  # DESeq2 
  if (opt$method == 'deseq2' | opt$method == 'all') {
    output_prefix <- change_dirs(opt$result_dir, 'DESeq2', comparison)
    relevel(dds$Condition, ref= combs[[c]][1])
    dds <- DESeq(dds)
    
    # Write filtered count matrix to file
    write.table(counts(dds), file=paste(output_prefix, "count_mtx_filtered.csv", sep=""), sep=",", quote=F, col.names=NA)
    
    # Write normalized values to file
    normalized_count_mtx_dds <- counts(dds, normalized=TRUE)
    write.table(normalized_count_mtx_dds, file=paste(output_prefix, "count_mtx_normalized.csv", sep=""), sep=",", quote=F, col.names=NA)
    
    # Obtain DEG results and output to file
    res_dds <- results(dds, contrast=c("Condition", combs[[c]][2], combs[[c]][1]))
    res_dds[['FoldChange']] <- 2^abs(res_dds[['log2FoldChange']])*(res_dds[['log2FoldChange']]/abs(res_dds[['log2FoldChange']]))
    write.table(res_dds[order(res_dds$log2FoldChange, decreasing=TRUE), ], file=paste(output_prefix, "DESeq2_FullResult_", comparison, ".csv", sep=""), sep=",", quote=F, col.names=NA)
    
  }
  
  # edgeR
  if (opt$method == 'edger' | opt$method == 'all') {
    output_prefix <- change_dirs(opt$result_dir, 'edgeR', comparison)
    write.table(edge$counts, file=paste(output_prefix, "count_mtx_filtered.csv", sep=""), sep=",", quote=F, col.names=NA)
    
    edge <- calcNormFactors(edge, method="TMM")
    normalized_count_mtx_edge <- cpm(edge, normalized.lib.sizes=TRUE) # method called cpm, but uses norm factors from above with TMM
    write.table(normalized_count_mtx_edge, file=paste(output_prefix, "count_mtx_normalized.csv", sep=""), sep=",", quote=F, col.names=NA)
  }
  
  # Filter out normalized low-expressed genes (bottom 10th percentile) by baseMean  REMOVE ALL < 10 baseMean
  if (opt$filter){
    
    # DESeq2
    if (opt$method == 'deseq2' | opt$method == 'all') {
      # Remove any genes below minimum baseMean
      output_prefix <- change_dirs(opt$result_dir, 'DESeq2', comparison)
      cat("DESeq2 genes and samples before filtering normalized genes:\n", dim(res_dds), "\n")
      cat("DESeq2 removing genes with baseMean expression <", opt$min_basemean, "\n")
      res_dds <- res_dds[res_dds$baseMean >= opt$min_basemean, ]
      cat("DESeq2 genes and samples after filtering normalized genes:\n", dim(res_dds), "\n")
    }
    
    # edgeR
    if (opt$method == 'edger' | opt$method == 'all') {
      output_prefix <- change_dirs(opt$result_dir, 'edgeR', comparison)
      cat("EdgeR genes and samples before filtering normalized genes:\n", dim(normalized_count_mtx_edge), "\n")
      cat("EdgeR removing genes with sum expression <", opt$min_basemean, "\n")
      norm_filtered_count_mtx_edge <- subset(normalized_count_mtx_edge, rowSums(normalized_count_mtx_edge) > opt$min_basemean)
      cat("EdgeR genes and samples after filtering normalized genes:\n", dim(norm_filtered_count_mtx_edge), "\n")
      
      # Write to file
      write.table(norm_filtered_count_mtx_edge, file=paste(output_prefix, "count_mtx_normalized_filtered.csv", sep=""), sep=",", quote=F, col.names=NA)
      
      # Update edgeR object with filtered genes
      edge <- DGEList(counts = norm_filtered_count_mtx_edge, group=mtx_info$Condition, remove.zeros=FALSE, norm.factors = edge$samples$norm.factors)
    }
  }
  
  # DESeq2
  if (opt$method == 'deseq2' | opt$method == 'all') {
    output_prefix <- change_dirs(opt$result_dir, 'DESeq2', comparison)
    
    # Generate plots
    png(paste(output_prefix, 'MAplot.png', sep=''))
    DESeq2::plotMA(res_dds, ylim=c(-2,2), ylab='log2FoldChange', alpha=0.1) 
    title(main = comparison, sub="blue if adjusted pvalue < 0.1", cex.sub=0.8)
    abline(h=c((0-opt$lfc), opt$lfc), col="red", lwd=2)
    invisible(capture.output(dev.off()))
    #invisible(capture.output(ggsave(paste(output_prefix, 'MAplot.png', sep=''), plot=plt)))
    
    resLFC <- lfcShrink(dds, coef=resultsNames(dds)[-1], type="apeglm")
    png(paste(output_prefix, 'shrunkMAplot.png', sep=''))
    DESeq2::plotMA(resLFC, ylim=c(-2,2), ylab='log2FoldChange', alpha=0.1)
    title(main = comparison, sub="blue if adjusted pvalue < 0.1", cex.sub=0.8)
    abline(h=c((0-opt$lfc), opt$lfc), col="red", lwd=2)
    invisible(capture.output(dev.off()))
    
    plt <- make_volcanoplot(res_dds, combs[[c]][1], combs[[c]][2], opt$lfc, 0.1)
    invisible(capture.output(ggsave(filename=paste(output_prefix, 'Volcanoplot.png', sep=''), plot=plt, dpi=320)))
    
    vsd <- vst(dds, blind=FALSE)
    #png(paste(output_prefix, 'PCAplot.png', sep=''))
    plt <- plotPCA(vsd, intgroup=c("Condition")) + 
      geom_text(aes(label = rownames(vsd@colData)), nudge_y = 0.2, size = 3)
    invisible(capture.output(ggsave(paste(output_prefix, 'PCAplot.png', sep=''), plot=plt)))
    #invisible(capture.output(dev.off()))
  }
  
  # edgeR
  if (opt$method == 'edger' | opt$method == 'all') {
    output_prefix <- change_dirs(opt$result_dir, 'edgeR', comparison)
    disp <- estimateDisp(edge)
    test <- exactTest(disp)
    res_edge <- topTags(test, n=dim(test$table)[1], adjust.method="BH", sort.by="logFC", p.value=1)
    
    # General plots
    png(paste(output_prefix, 'MDSplot.png', sep=''))
    plotMDS(edge, col=as.numeric(edge$samples$group), plot=TRUE)
    invisible(capture.output(dev.off()))
    
    png(paste(output_prefix, 'BCVplot.png', sep=''))
    plotBCV(disp)
    invisible(capture.output(dev.off()))
    
    de.tags <- subset(res_edge$table, res_edge$table$FDR < 0.1)
    png(paste(output_prefix, 'Smearplot.png', sep=''))
    plotSmear(test, ylim=c(-2,2), ylab='log2foldchange', de.tags=rownames(de.tags), n=dim(test)[1])
    title(main = comparison, sub="red if adjusted pvalue < 0.1", cex.sub=0.8)
    abline(h=c((0-opt$lfc), opt$lfc), col="blue", lwd=2)
    invisible(capture.output(dev.off()))
    
    colnames(res_edge$table) <- c('log2FoldChange', 'log2CPM', 'pvalue', 'padj')
    res_edge$table[['FoldChange']] <- 2^abs(res_edge$table[['log2FoldChange']])*(res_edge$table[['log2FoldChange']]/abs(res_edge$table[['log2FoldChange']]))
    res_edge <- as.data.frame(res_edge)
    
    write.table(res_edge[order(res_edge$log2FoldChange, decreasing=TRUE), ], file=paste(output_prefix, "edgeR_FullResult_", comparison, ".csv", sep=""), sep=",", quote=F, col.names=NA)
  }
  
  # Analyze results from methods
  results <- list()
  if (opt$method == 'deseq2' | opt$method == 'all') {
    results[['DESeq2']] <- res_dds
  }
  if (opt$method == 'edger' | opt$method == 'all') {
    results[['edgeR']] <- res_edge
  }
  
  for (r in names(results)){
    
    if (opt$method == tolower(r) | opt$method == 'all') {
      
      cat("\nAnalyzing", r, "results...\n")
    
      output_prefix <- change_dirs(opt$result_dir, r, comparison)
      res <- results[[r]]
    
      # ========== Organize DEGs ==============
      out_dirs <- list()
      out_dirs[combs[[c]][1]] <- paste(output_prefix, combs[[c]][1], '/', sep='')
      out_dirs[combs[[c]][2]] <- paste(output_prefix, combs[[c]][2], '/', sep='')
      out_dirs["DEGs"] <- paste(output_prefix, "DEGs", '/', sep='')
      if (!file.exists(out_dirs[[1]])) {
        dir.create(out_dirs[[1]])
      }
      if (!file.exists(out_dirs[[2]])) {
        dir.create(out_dirs[[2]])
      }
      if (!file.exists(out_dirs[[3]])) {
        dir.create(out_dirs[[3]])
      }
      
      # Remove any NAs
      res <- res[rowSums(is.na(res)) == 0, ]
      
      # Get significantly Upregulated (log2FoldChange > 0)
      res_up <- subset(res, log2FoldChange >= opt$lfc & pvalue <= opt$pvalue)
      res_up_sorted <- res_up[order(res_up$log2FoldChange, decreasing=TRUE),]
      write.table(res_up_sorted, file=paste(out_dirs[[1]], r, "_Result_", combs[[c]][1], ".csv", sep=""), sep=",", quote=F, col.names=NA)
      
      # Get significantly Downregulated (log2FoldChange < 0)
      res_down <- subset(res, log2FoldChange <= (0 - opt$lfc) & pvalue <= opt$pvalue)
      res_down_sorted <- res_down[order(res_down$log2FoldChange, decreasing=TRUE),]
      write.table(res_down_sorted, file=paste(out_dirs[[2]], r, "_Result_", combs[[c]][2], ".csv", sep=""), sep=",", quote=F, col.names=NA)
      
      # Get all significantly changed genes
      res_changed <- subset(res, (log2FoldChange >= opt$lfc | log2FoldChange <= (0 - opt$lfc)) & pvalue <= opt$pvalue)
      res_changed_sorted <- res_changed[order(abs(res_changed$log2FoldChange), decreasing=TRUE),]
      write.table(res_changed_sorted, file=paste(out_dirs[["DEGs"]], r, "_Result_DEGs", ".csv", sep=""), sep=",", quote=F, col.names=NA)
      
      # Heatmap
      plt <- make_heatmapplot(res_changed_sorted, combs[[c]][1], combs[[c]][2], n=40)
      invisible(capture.output(ggsave(filename=paste(out_dirs[["DEGs"]], 'Heatmapplot.png', sep=''), plot=plt, dpi=320)))
      
      genes <- list()
      if (dim(res_up_sorted)[1] != 0){
        genes[[combs[[c]][1]]] <- rownames(res_up_sorted)
      }
      if (dim(res_down_sorted)[1] != 0){
        genes[[combs[[c]][2]]] <- rownames(res_down_sorted)
      }
      if (dim(res_changed_sorted)[1] != 0){
        genes[["DEGs"]] <- rownames(res_changed_sorted)
      }
      
      genes_entrez <- list()
      for (n in names(genes)){
        if (opt$assembly == "hg19" | opt$assembly == "hg38"){
          genes_entrez[[n]] <- mapIds(org.Hs.eg.db, keys = genes[[n]], column = "ENTREZID", keytype = "SYMBOL")
        }else if (opt$assembly == "mm9" | opt$assembly == "mm10"){
          genes_entrez[[n]] <- mapIds(org.Mm.eg.db, keys = genes[[n]], column = "ENTREZID", keytype = "SYMBOL")
        }else if (opt$assembly == "rn6"){
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
                                       OrgDb=anno_ref$annoDb,
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
                invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'GO_', ont, '_', n, '.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                # Write annotations to csv
                write.table(as.data.frame(compGO), file=paste(out_dirs[[n]], 'GO_', ont, '_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_pheatmapplot(compGO@compareClusterResult, res, assembly=opt$assembly, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE, title=paste("GO (", ont, ") - ", n, sep=""), ylabel="GO Term")
                invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'GO_', ont, '_', n, '_pheatmap_by_gene.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                plt <- make_pheatmapplot(compGO@compareClusterResult, res, assembly=opt$assembly, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE, title=paste("GO (", ont, ") - ", n, sep=""), ylabel="GO Term")
                invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'GO_', ont, '_', n, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
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
                                       organism=anno_ref$keggOrg
            ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
            if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
              #plt <- dotplot(compKEGG, showCategory = 8, title = paste("KEGG -", n, sep=""))
              plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', n, sep=""), ylabel="KEGG Category", colour=colour, n=15)
              invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'KEGG_annotation_', n, '.png', sep=''), plot=plt, dpi=320)))
              remove(plt)
              # Write annotations to csv
              write.table(as.data.frame(compKEGG), file=paste(out_dirs[[n]], 'KEGG_annotation_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', n, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE, ylabel="KEGG Category")
              invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'KEGG_annotation_', n, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
              remove(plt)
              plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', n, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE, ylabel="KEGG Category")
              invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'KEGG_annotation_', n, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
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
                            OrgDb = anno_ref$annoDb, 
                            pAdjustMethod = "BH"
              )
              if ((!is.null(gsea)) & (dim(gsea@result)[1] > 0)){
                gsea@result <- gsea@result[order(gsea@result$p.adjust, decreasing=FALSE),] # Sort by most signiicant
                gsea@result <- head(gsea@result, n=20) # retain only top to plot
                plt <- dotplot(gsea, showCategory = 20, title = paste("GSEA (", ont, ") - ", n, sep=""))
                df <- plt$data
                plt <- make_dotplot(df, title=paste("GSEA (", ont, ") - ", n, sep=""), ylabel="GSEA", colour=colour, n=15)
                invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'GSEA_', ont, '_', n, '.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                
                # Write annotations to csv
                write.table(as.data.frame(gsea$result), file=paste(out_dirs[[n]], 'GSEA_', ont, '_', n, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_pheatmapplot(df, res, assembly=opt$assembly, title=paste("GSEA (", ont, ") - ", n, sep=""), ylabel="GSEA", heat_colour=heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE)
                invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'GSEA_', ont, '_', n, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                plt <- make_pheatmapplot(df, res, assembly=opt$assembly, title=paste("GSEA (", ont, ") - ", n, sep=""), ylabel="GSEA", heat_colour=heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE)
                invisible(capture.output(ggsave(filename=paste(out_dirs[[n]], 'GSEA_', ont, '_', n, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
              }
              
            }
          },error = function(e)
          {
            message(e)
          }
        ) 
      }
    
    }else{
      cat("\nMethod not supported:", r)
    }
    cat("Completed", r,  "\n")
  }
  
}
