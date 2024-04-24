# -----------------------------------------------------------
# Script Name: Analyze Peaks
# Purpose: This script performs differential peaks analysis using DiffBind.
# Author: Eric Arezza
# Date: 2024-04-23
# Version: 1.0
# -----------------------------------------------------------

# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
  library(RDAVIDWebService)
  library(DiffBind)
  library(tidyverse)
  library(ChIPseeker)
  library(ReactomePA)
  library(clusterProfiler)
  library(vulcan)
  library(UpSetR)
  library(dplyr)
  library(VennDiagram)
  library(eulerr)
  require(gridExtra)
  library(RColorBrewer)
  library(pheatmap)
  library(grid)
  library(optparse)
}))

# ======= Get command-line optional arguments =======
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="DiffBind-formatted sample sheet", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default="mm10", help="Assembly to annotate genes/peaks (e.g. hg19, hg38, mm9, mm10, rn6)", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="Peaks_Analysis/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-d", "--database"), type="character", default="ucsc", help="Database reference for peaks gene annotations, ucsc (default) or ensembl", metavar="character"),
  make_option(c("-l", "--annotation_level"), type="character", default="transcript", help="Level parameter for annotatePeak, 'gene' or 'transcript'", metavar="character"),
  make_option(c("--combine_callers"), type="logical", action="store_true", default=FALSE, help="Flag to add peaks from callers instead of taking consensus peaks", metavar="logical"),
  make_option(c("--combine_replicates"), type="logical", action="store_true", default=FALSE, help="Flag to add peaks from replicates instead of taking consensus peaks", metavar="logical"),
  make_option(c("-b", "--blacklisted_keep"), type="logical", action="store_true", default=FALSE, help="Flag to keep blacklisted regions in raw peaks files", metavar="logical"),
  make_option(c("--lfc"), type="double", default=0.585, help="Magnitude of log2foldchange to define significant up/down regulation of genes", metavar="double"),
  make_option(c("--fdr"), type="double", default=0.05, help="Significance threshold (false discovery rate, a.k.a. p.adjust value) for DEGs", metavar="double"),
  make_option(c("--occupancy_only"), type="logical", action="store_true", default=FALSE, help="Flag to only perform peaks occupancy analysis", metavar="logical"),
  make_option(c("--david_user"), type="character", default="earezza@ohri.ca", help="User email for DAVID web tools (must be registered, https://david.ncifcrf.gov/content.jsp?file=DAVID_WebService.html)", metavar="character"),
  make_option(c("--minGSSize"), type="integer", default=10, help="minimal size of genes annotated for testing", metavar="integer"),
  make_option(c("--maxGSSize"), type="integer", default=500, help="maximal size of genes annotated for testing", metavar="integer")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cat("Command-line options:\n")
for (i in which(names(opt) != "help")) {
  cat(names(opt)[i], '=', paste(opt)[i], "\n")
}
cat("log2FC of", opt$lfc, "equates to FC of", round(2^0.585, 2), '\n')

if (str_sub(opt$result_dir, -1) != "/"){
  opt$result_dir = cat(opt$result_dir, '/', sep='')
}
if (!file.exists(opt$file)){
  cat(opt$file, "does not exist...check path and current working directory.")
  q()
}
if (!(opt$assembly %in% c('mm10', 'mm9', 'hg38', 'hg19', 'rn6'))){
  cat(opt$assembly, "not a valid choice. Only supports mm9, mm10, hg19, hg38, rn6 assemblies.")
  q()
}


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
    while (i < length(strsplit(s, ' ')[[1]]) + 1){
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
    theme_classic() +
    theme(axis.text.y = element_text(colour=rev(head(df$ycolour, n=n)))) +
    scale_color_gradient(low = "black", high = colour) +
    ggtitle(title) 
  plt$labels$x <- "-log(p.adjust)"
  plt$labels$y <- ylabel
  plt$labels$size <- "GenePercentage"
  plt$labels$colour <- "GeneCount"
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
      }else if (assembly == "mm10" | assembly == "mm9"){
        gene_group <- mapIds(org.Mm.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }else if (assembly == "rn6"){
        gene_group <- mapIds(org.Rn.eg.db, keys = gene_group, column = "SYMBOL", keytype = "ENTREZID")
      }
    }
    # For RNAseq degs
    #d[gene_group, a] <- res[gene_group, ]$log2FoldChange
    
    # For Peaks degs
    temp <- res[res$SYMBOL %in% gene_group, ]
    temp <- aggregate(temp, list(temp$SYMBOL), mean)
    row.names(temp) <- temp$Group.1
    d[temp$Group.1, a] <- temp$log2FoldChange
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

# # ========= START OCCUPANCY ANALYSIS =========
# Here, peaks declared by peak caller(s) are used to identify
# differential expression between conditions. Overlapping peaks 
# between replicates/conditions are determined by the range of the peaks. 
# In this case, using occupancy alone for DE provides a less conservative
# analysis for DE by simply considering where peaks exist.
#
# Following this, affinity analysis below can provide a more conservative
# result for DE since the read counts are accounted for (peak shapes).


# ========= SETUP RUN AND VARIABLES =========
# ========= Load peaksets =========
setwd(dirname(opt$file))
samplesheet <- basename(opt$file)

# Set directories
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}
cat("Output files will be in", opt$result_dir, "\n")

result_dir <- paste(opt$result_dir, "Occupancy_Analysis/", sep='')
if (!file.exists(result_dir)) {
  dir.create(result_dir)
}

supplementary_dir <- paste(result_dir, "Supplementary/", sep='')
if (!file.exists(supplementary_dir)) {
  dir.create(supplementary_dir)
}

#output_prefix <- gsub('.csv', '_', paste(result_dir, samplesheet, sep=""))
#output_prefix <- gsub("diffbind_samplesheet_", "", output_prefix)

# Get annotations reference and respective promoter regions
anno_ref <- load_annotation(opt$assembly, opt$database)
promoters <- getPromoters(TxDb=anno_ref$txdb, upstream=3000, downstream=3000)

# Average for each sample
fragment_size <- 1:length(read.csv(samplesheet)$SampleID)
frag_sizes <- list()
for (b in unique(read.csv(samplesheet)$Condition)){
  for (r in unique(read.csv(samplesheet)$Replicate)){
    png(paste(supplementary_dir, 'fragment_length_', b, '-', r, '.png', sep=""))
    mean_fragment_size <- average_fragment_length(read.csv(samplesheet)$bamReads[[which(read.csv(samplesheet)$Condition == b)[1]]], plot=TRUE)
    for (i in which(read.csv(samplesheet)$Condition == b & read.csv(samplesheet)$Replicate == r)){
      fragment_size[i] <- mean_fragment_size
      frag_sizes[[b]] <- mean_fragment_size
    }
    invisible(capture.output(dev.off())) 
  }
}
invisible(capture.output(gc())) 
#fragment_size <- 125 # default

dbObj <- dba(sampleSheet=samplesheet, minOverlap=1,
             config=data.frame(th=opt$fdr,
                               #DataType=DBA_DATA_GRANGES, 
                               RunParallel=TRUE,
                               minQCth=15, 
                               fragmentSize=fragment_size,
                               reportInit="DBA",
                               bUsePval=FALSE
             )
)
cat("Raw peaksets:\n")
dbObj

# Colour codes for consistency in plots (add more colours if needed)
colours <- c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF", "#e69e02", "#00A9FF", "#C77CFF", "#FF61CC", 
             "#FF0000", "#FF4D00", "#80FF00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2" , "#D55E00" , "#CC79A7"
             )
conditions_colour_code <- list()
for (i in 1:length(unique(dbObj$samples$Condition))) {
  conditions_colour_code[[unique(dbObj$samples$Condition)[i]]] <- colours[i]
}
conditions_colour_code[['Shared']] <- "grey"

png(filename=paste(supplementary_dir, 'raw_heatmap.png', sep=''))
dba.plotHeatmap(dbObj)
invisible(capture.output(dev.off()))
invisible(capture.output(gc()))

# Show overlap rates for each condition
cat("Peak overlaps in at least (1, 2, ...) replicates/callers for each condition:\n")
png(paste(supplementary_dir, "raw_overlap_rates.png", sep=""))
par(mfrow=c(length(unique(dbObj$samples$Condition)), 1), mar = c(2, 4, 4, 2))
for (c in unique(dbObj$samples$Condition)) {
  cat('\n', c, '\n')
  olap.rate <- dba.overlap(dbObj, dba.mask(dbObj, attribute=DBA_CONDITION, value=c, combine='or'), mode=DBA_OLAP_RATE)
  cat(olap.rate)
  cat('\n')
  plot(olap.rate, type='l', xlim=c(0.75, length(olap.rate)), ylab='# overlapping peaks', xlab='# peaksets (replicates and peak callers)', col=conditions_colour_code[[c]])
  text(x=1:length(olap.rate), y=olap.rate, olap.rate, cex=0.75)
  axis(side=1, at=1:length(olap.rate))
  title(main=c("Overlap Rate", c))
}
invisible(capture.output(dev.off()))
invisible(capture.output(gc()))

# ========= Remove Blacklisted Regions =========
# Remove blacklisted regions to ignore irrelevant peaks (blacklisted regions from ENCODE, genome selected is based on prediction from bam files)
tryCatch (
  {
    if (!opt$blacklisted_keep){
      dbObj.noblacklist <- dba.blacklist(dbObj, blacklist=TRUE, greylist=FALSE)
      blacklisted_peaks <- dba.blacklist(dbObj.noblacklist, Retrieve=DBA_BLACKLISTED_PEAKS)
      cat("After blacklisted regions removed:\n")
      dbObj.noblacklist
    }else{
      cat("\nBlacklisted regions not removed, proceeding with raw peaksets...\n")
      dbObj.noblacklist <- dbObj
    }
  },error = function(e)
  {
    message(e)
  }
)
if (!exists("dbObj.noblacklist")) {
  cat("\nBlacklisted regions not removed, proceeding with raw peaksets...\n")
  dbObj.noblacklist <- dbObj
  dbObj.noblacklist
}

png(paste(supplementary_dir, 'raw_noblacklist_heatmap.png', sep=''))
dba.plotHeatmap(dbObj.noblacklist)
invisible(capture.output( dev.off() ))
invisible(capture.output(gc()))

# Show overlap rates for each condition
cat("Peak overlaps in at least (1, 2, ...) replicates/callers for each condition:\n")
png(paste(supplementary_dir, "raw_noblacklist_overlap_rates.png", sep=""))
par(mfrow=c(length(unique(dbObj$samples$Condition)), 1), mar = c(2, 4, 4, 2))
for (c in unique(dbObj$samples$Condition)) {
  cat('\n', c, '\n')
  olap.rate <- dba.overlap(dbObj, dba.mask(dbObj.noblacklist, attribute=DBA_CONDITION, value=c, combine='or'), mode=DBA_OLAP_RATE)
  cat(olap.rate)
  cat('\n')
  plot(olap.rate, type='l', xlim=c(0.75, length(olap.rate)), ylab='# overlapping peaks', xlab='# peaksets (replicates and peak callers)', col=conditions_colour_code[[c]])
  text(x=1:length(olap.rate), y=olap.rate, olap.rate, cex=0.75)
  axis(side=1, at=1:length(olap.rate))
  title(main=c("Overlap Rate", c))
}
invisible(capture.output(dev.off()))
invisible(capture.output(gc()))

# ========= Get Consensus Peaks =========
if (opt$combine_replicates == TRUE){
  rep_overlaps <- 1 # add all peaks from replicates
}else{
  rep_overlaps <- 2 # add only consensus peaks (peaks must be in at least 2 replicates)
}

# If using peaks from multiple peak callers (defined in Factor column of samplesheet)
if (length(unique(dbObj$samples$Factor)) > 1){
  # Get consensus overlaps in peaksets for each condition (peaks must be overlapping in majority (2/3rds) of peak callers) or add all callers' peaks
  if (opt$combine_callers == TRUE){
    dbObj.total <- dba.peakset(dbObj.noblacklist, consensus=c(DBA_CONDITION, DBA_REPLICATE), minOverlap=1)
  }else{
    dbObj.total <- dba.peakset(dbObj.noblacklist, consensus=c(DBA_CONDITION, DBA_REPLICATE), minOverlap=0.66)
  }
  # resulting consensus between callers
  dbObj.caller_consensus <- dba(dbObj.total, mask=dbObj.total$masks$Consensus, minOverlap=1)
  if (length(unique(dbObj$samples$Replicate)) > 1){
    # Get consensus between replicates for each condition (peaks must be overlapping in at least 2 replicates)
    dbObj.final <- dba.peakset(dbObj.caller_consensus, consensus=c(DBA_CONDITION), minOverlap=rep_overlaps)
    maskname <- names(dbObj.final$masks)[grepl('Replicate.1-2', names(dbObj.final$masks))]
    # resulting consensus between replicates
    if (length(maskname) == 1){
      dbObj.consensus <- dba(dbObj.final, mask=dbObj.final$masks[[maskname]], minOverlap=1)
    }else if (length(maskname) == 2){
      dbObj.consensus <- dba(dbObj.final, mask=(dbObj.final$masks[[maskname[1]]] | dbObj.final$masks[[maskname[2]]]), minOverlap=1)
    }
  }else{
    dbObj.consensus <- dbObj.caller_consensus
  }
}else{
  if (length(unique(dbObj$samples$Replicate)) > 1){
    # consensus between replicates
    dbObj.total <- dba.peakset(dbObj.noblacklist, consensus=c(DBA_CONDITION), minOverlap=rep_overlaps)
    maskname <- names(dbObj.total$masks)[grepl('Replicate.1-2', names(dbObj.total$masks))]
    if (length(maskname) == 1){
      dbObj.consensus <- dba(dbObj.total, mask=dbObj.total$masks[[maskname]], minOverlap=1)
      dbObj.caller_consensus <- dbObj.consensus
    }else if (length(maskname) == 2){
      dbObj.consensus <- dba(dbObj.total, mask=(dbObj.total$masks[[maskname[1]]] | dbObj.total$masks[[maskname[2]]]), minOverlap=1)
      dbObj.caller_consensus <- dbObj.consensus
    }
  } else {
    # When only 1 peak caller and 1 replicate available
    dbObj.total <- dbObj.noblacklist
    dbObj.consensus <- dbObj.noblacklist
    dbObj.caller_consensus <- dbObj.noblacklist
    
  }
}
dbObj.caller_consensus
dbObj.consensus


# Add fragment sizes to new objects (helps to later run affinity analysis)
for (c in unique(dbObj$samples$Condition)){
  for (i in 1:length(dbObj.caller_consensus$mask[[c]])){
    if (dbObj.caller_consensus$mask[[c]][[i]]){
      dbObj.caller_consensus$config$fragmentSize[i] <- frag_sizes[[c]]
    }
  }
}
for (c in unique(dbObj$samples$Condition)){
  for (i in 1:length(dbObj.consensus$mask[[c]])){
    if (dbObj.consensus$mask[[c]][[i]]){
      dbObj.consensus$config$fragmentSize[i] <- frag_sizes[[c]]
    }
  }
}

# Re-sort colours if condition orders changed after consensus (occurs when one condition has only 1 replicate, consensus must be added manually...)
# dbObj.consensus <- dba(dbObj.final, mask=(dbObj.final$masks$`Replicate.1-2` | dbObj.final$masks$CONDITION_WITH_ONE_REPLICATE), minOverlap=1)
temp <- list()
for (i in 1:length(conditions_colour_code)){
  temp[names(conditions_colour_code[which(names(conditions_colour_code) == dba.show(dbObj.consensus)$Condition[i])])] <- conditions_colour_code[which(names(conditions_colour_code) == dba.show(dbObj.consensus)$Condition[i])]
}
conditions_colour_code <- temp

# Consensus peaks from all conditions (all relevant peaks)
consensus_peaks <- dba.peakset(dbObj.consensus, bRetrieve=TRUE)

result_dirs <- list()
for (p in unique(dba.show(dbObj.consensus)$Condition)){
  result_dirs[[p]] <- paste(result_dir, p, "/", sep='')
  if (!file.exists(result_dirs[[p]])) {
    dir.create(result_dirs[[p]])
  }
}

# Output consensus peaksets
raw_peaks <- list()
for (c in unique(dba.show(dbObj.consensus)$Condition)){
  p <- dba.peakset(dbObj.consensus, dbObj.consensus$masks[[c]], bRetrieve=TRUE)
  write.table(as.data.frame(p), file=paste(result_dirs[[c]], c, '_consensus.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
  raw_peaks[[c]] <- p
}

# ========= Get Annotations =========
peakAnnoList <- list()
for (p in names(raw_peaks)){
  cat("\nAnnotating", p, ' consensus\n')
  anno <- annotatePeak(raw_peaks[[p]], 
                       TxDb=anno_ref$txdb,
                       annoDb=anno_ref$annoDb,
                       level=opt$annotation_level,
                       tssRegion=c(-3000, 3000)
  )
  peakAnnoList[[p]] <- anno
  write.table(anno@anno, file=paste(result_dirs[[p]], c, '_consensus_annotated.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
}

# Plot UpSet
upsetlist <- list()
for (p in names(peakAnnoList)) {
  upsetlist[[p]] <- peakAnnoList[[p]]@anno$SYMBOL
}
upset_colors <- list()
for (n in names(peakAnnoList)){
  upset_colors[[ conditions_colour_code[[n]] ]] <- length(unique(peakAnnoList[[n]]@anno$SYMBOL))
}
upset_colors <- sort(unlist(upset_colors), decreasing=TRUE)
png(paste(result_dir, 'consensus_annotated-genes_upsetplot.png', sep=''),
    width = 1920,
    height = 1080,
    res=200
    )
upset(fromList(upsetlist), 
             order.by = "freq", 
             nsets = length(names(peakAnnoList)),
      sets.bar.color = names(upset_colors),
      empty.intersections = "on",
      set_size.show = TRUE,
      set_size.angles = 0,
      set_size.scale_max = dim(fromList(upsetlist))[[1]],
      sets.x.label = "Gene Set Size",
      mainbar.y.label = "Intersection Size of Gene Sets",
      mb.ratio = c(0.7, 0.3)
      ) 
grid.text("Consensus Peaksets Annotated Genes",x = 0.65, y=0.95, gp=gpar(fontsize=10))
dev.off()
invisible(capture.output(gc()))

          
# Differentially bound peaks between conditions (peaks unique to each condition)
differential_peaks <- dba.overlap(dbObj.consensus, dbObj.consensus$masks$All, DataType=DBA_DATA_GRANGES)
unique_peaks <- list()
for (c in unique(dba.show(dbObj.consensus)$Condition)){
  unique_peaks[[c]] <- differential_peaks[[which(unique(dba.show(dbObj.consensus)$Condition) == c)]]
  cat("\n", c, "has", length(unique_peaks[[c]]), "unique peaks.\n")
}
# Peaks shared between all conditions
shared_peaks <- list()
shared_peaks[["Shared"]] <- differential_peaks$inAll
cat("\n", length(shared_peaks[["Shared"]]), "all shared peaks.\n")

if (length(unique(dbObj$samples$Condition)) == 3){
  for (c in unique(dba.show(dbObj.consensus)$Condition)){
    i <- which(unique(dba.show(dbObj.consensus)$Condition) == c) + 3
    pair <- unique(dba.show(dbObj.consensus)$Condition)[which(unique(dba.show(dbObj.consensus)$Condition) != c)]
    unique_peaks[[paste(pair, collapse='_and_')]] <- differential_peaks[[which(unique(dba.show(dbObj.consensus)$Condition) == c) + 3]]
    conditions_colour_code[[paste(pair, collapse='_and_')]] <- colours[i]
    cat("\n", paste(pair, collapse='_and_'), "have", length(unique_peaks[[paste(pair, collapse='_and_')]]), "shared peaks.\n")
  }
}

if (length(unique(dbObj$samples$Condition)) == 4){
  e = c(
    "A"=length(unique_peaks[[names(unique_peaks)[1]]]), 
    "B"=length(unique_peaks[[names(unique_peaks)[2]]]),
    "C"=length(unique_peaks[[names(unique_peaks)[3]]]),
    "D"=length(unique_peaks[[names(unique_peaks)[4]]]),
    "A&B"=length(differential_peaks$AandB),
    "A&C"=length(differential_peaks$AandC),
    "A&D"=length(differential_peaks$AandD),
    "B&C"=length(differential_peaks$BandC),
    "B&D"=length(differential_peaks$BandD),
    "C&D"=length(differential_peaks$CandD),
    "A&B&C"=length(differential_peaks$notD),
    "A&B&D"=length(differential_peaks$notC),
    "A&C&D"=length(differential_peaks$notB),
    "B&C&D"=length(differential_peaks$notA),
    "A&B&C&D"=length(shared_peaks[[names(shared_peaks)[1]]])
  )
  names(e) = c(names(unique_peaks)[1], names(unique_peaks)[2], names(unique_peaks)[3], names(unique_peaks)[4],
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[3] , '&' , names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&', names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&', names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[3] , '&', names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[3] , '&', names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&', names(unique_peaks)[3] , '&', names(unique_peaks)[4], sep='')
  )
  for (c in (length(unique(dbObj$samples$Condition))+1):(length(names(e))-1) ){
    combo <- str_replace_all(names(e)[c], "&", "_and_")
    unique_peaks[[combo]] <- differential_peaks[[c, ]]
    conditions_colour_code[[combo]] <- colours[i]
    cat("\n", combo, "have", length(unique_peaks[[combo]]), "shared peaks.\n")
  }
}
conditions_colour_code[['Shared']] <- "grey"



# Plots
if (length(unique(dba.show(dbObj.consensus)$Condition)) == 2){
  grid.newpage()
  g = draw.pairwise.venn(area1=length(unique_peaks[[names(unique_peaks)[1]]])+length(shared_peaks[[names(shared_peaks)[1]]]), 
                         area2=length(unique_peaks[[names(unique_peaks)[2]]])+length(shared_peaks[[names(shared_peaks)[1]]]),
                         cross.area=length(shared_peaks[[names(shared_peaks)[1]]]),
                         category=names(unique_peaks),
                         fill=unname(unlist(conditions_colour_code))[1:2],
                         col=NA,
                         #cat.pos=c(0,0),
                         cat.dist = c(0,0))
  plt <- grid.arrange(gTree(children=g), top="Binding Site Overlaps", bottom=gsub('/', '', opt$result_dir))
  invisible(capture.output(ggsave(filename=paste(result_dir, 'consensus_peaks_venn.png', sep=''), plot=plt)))
  rm(g)
}else if (length(unique(dba.show(dbObj.consensus)$Condition)) == 3){
  e = c(
    "A"=length(unique_peaks[[names(unique_peaks)[1]]]), 
    "B"=length(unique_peaks[[names(unique_peaks)[2]]]),
    "C"=length(unique_peaks[[names(unique_peaks)[3]]]),
    "A&B"=length(differential_peaks$notC),
    "B&C"=length(differential_peaks$notA),
    "A&C"=length(differential_peaks$notB),
    "A&B&C"=length(shared_peaks[[names(shared_peaks)[1]]])
  )
  names(e) = c(names(unique_peaks)[1], names(unique_peaks)[2], names(unique_peaks)[3],
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&' , names(unique_peaks)[3], sep='')
  )
  
  #png(paste(output_prefix, 'raw_consensus_peaks.png', sep=""))
  plt <- plot(euler(e), main=gsub('/', '', opt$result_dir), quantities=TRUE, fills=unname(unlist(conditions_colour_code)))
  invisible(capture.output(ggsave(filename=paste(result_dir, 'consensus_peaks_venn.png', sep=''), plot=plt)))
  #invisible(capture.output(dev.off()))
  rm(e)
}else if (length(unique(dba.show(dbObj.consensus)$Condition)) == 4){
  e = c(
    "A"=length(unique_peaks[[names(unique_peaks)[1]]]), 
    "B"=length(unique_peaks[[names(unique_peaks)[2]]]),
    "C"=length(unique_peaks[[names(unique_peaks)[3]]]),
    "D"=length(unique_peaks[[names(unique_peaks)[4]]]),
    "A&B"=length(differential_peaks$AandB),
    "A&C"=length(differential_peaks$AandC),
    "A&D"=length(differential_peaks$AandD),
    "B&C"=length(differential_peaks$BandC),
    "B&D"=length(differential_peaks$BandD),
    "C&D"=length(differential_peaks$CandD),
    "A&B&C"=length(differential_peaks$notD),
    "A&B&D"=length(differential_peaks$notC),
    "A&C&D"=length(differential_peaks$notB),
    "B&C&D"=length(differential_peaks$notA),
    "A&B&C&D"=length(shared_peaks[[names(shared_peaks)[1]]])
  )
  names(e) = c(names(unique_peaks)[1], names(unique_peaks)[2], names(unique_peaks)[3], names(unique_peaks)[4],
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[3] , '&' , names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&', names(unique_peaks)[3], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&', names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[3] , '&', names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[2] , '&' , names(unique_peaks)[3] , '&', names(unique_peaks)[4], sep=''),
               paste(names(unique_peaks)[1] , '&' , names(unique_peaks)[2] , '&', names(unique_peaks)[3] , '&', names(unique_peaks)[4], sep='')
  )
  
  #png(paste(output_prefix, 'raw_consensus_peaks.png', sep=""))
  plt <- plot(euler(e), main=gsub('/', '', opt$result_dir), quantities=TRUE, fills=unname(unlist(conditions_colour_code)))
  invisible(capture.output(ggsave(filename=paste(result_dir, 'consensus_peaks_venn.png', sep=''), plot=plt)))
  #invisible(capture.output(dev.off()))
  rm(e)
}

invisible(capture.output(gc()))

png(paste(supplementary_dir, 'consensus_heatmap.png', sep=''))
dba.plotHeatmap(dbObj.consensus)
invisible(capture.output( dev.off() ))
invisible(capture.output(gc()))


tryCatch(
  {
    plt <- dba.plotPCA(dbObj, masks=!dbObj.total$masks$Consensus, attributes=DBA_CONDITION, label=DBA_ID, vColors=(colours))
    #plt$main <- "PCA"
    invisible(capture.output(ggsave(filename=paste(supplementary_dir, 'pca_condition.png', sep=''), plot=grid.arrange(plt))))
    
    if(length(unique(dbObj$samples$Factor)) > 1){
      plt <- dba.plotPCA(dbObj, masks=!dbObj.total$masks$Consensus, attributes=DBA_FACTOR, label=DBA_ID)
      invisible(capture.output( ggsave(filename=paste(supplementary_dir, 'pca_factor.png', sep=''), plot=grid.arrange(plt)) ))
      #invisible(capture.output(dev.off()))
    }
  },error = function(e)
  {
    message(e)
  }
)
invisible(capture.output(gc()))

# Plot peaks over genome
tryCatch(
  {
    plt <- covplot(c(unique_peaks[1:length(unique(dbObj$samples$Condition))], shared_peaks), title="Peaks over Genome") + 
      scale_color_manual(values=rev(c(unlist(unname(conditions_colour_code[1:length(unique(dbObj$samples$Condition))])), 'grey'))) + 
      scale_fill_manual(values=rev(c(unlist(unname(conditions_colour_code[1:length(unique(dbObj$samples$Condition))])), 'grey')))
    invisible(capture.output(ggsave(filename=paste(result_dir, 'genome_peaks.png', sep=''), plot=plt, dpi=320)))
    plt <- plt + facet_grid(chr ~ .id)
    invisible(capture.output(ggsave(filename=paste(result_dir, 'genome_peaks_split.png', sep=''), plot=plt, dpi=320)))
    rm(plt)
  },error = function(e)
  {
    message(e)
  }
)
invisible(capture.output(gc()))

peaks <- c(unique_peaks, shared_peaks)

for (p in names(peaks)){
  result_dirs[[p]] <- paste(result_dir, p, "/", sep='')
  if (!file.exists(result_dirs[[p]])) {
    dir.create(result_dirs[[p]])
  }
}

# Output peaks to bed files
for (p in names(peaks)){
  write.table(as.data.frame(peaks[[p]]), file=paste(result_dirs[[p]], p, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
}

# Plot peaks related to TSS sites
tryCatch(
  {
    tagMatrices <- list()
    for (p in names(peaks)){
      tagMatrix <- getTagMatrix(peaks[[p]], windows=promoters)
      if (length(tagMatrix) == 0){
        cat("No peaks at promoter sites for", p, "\n")
        rm(tagMatrix)
        #break
      }else{
        tagMatrices[[p]] <- tagMatrix
        cat(dim(tagMatrix)[[1]], "peaks at promoter sites for", p, "\n")
        plt <- tagHeatmap(tagMatrix, 
                          xlab="bp at TSS", 
                          ylab="Peaks", 
                          title=paste(dim(tagMatrix)[[1]],'Peaks at Promoters', p, sep=" - "),
                          palette=if_else(conditions_colour_code[[p]] == "#00BFC4", 'Greens', 'Reds'), 
        )
        invisible(capture.output(ggsave(paste(result_dirs[[p]], 'TSS_heatmap_', p, '_peaks.png', sep=''), plot=plt, dpi=320)))
        rm(tagMatrix)
        invisible(capture.output(gc()))
      }
    }
  },error = function(e)
  {
    message(e)
  }
)
invisible(capture.output(gc()))

# Plot TSS profile of peaks
tryCatch(
  {
    if (length(unique(dba.show(dbObj.consensus)$Condition)) == 3){
      plt <- plotAvgProf(tagMatrices[1:3], xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
        scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[1:3]) +
        scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[1:3])
      invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_unique-peaks.png', sep=''), plot=plt, dpi=320)))
      invisible(capture.output(gc()))
      plt <- plotAvgProf(tagMatrices[4:6], xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
        scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[4:6]) +
        scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[4:6])
      invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_pairs-peaks.png', sep=''), plot=plt, dpi=320)))
      invisible(capture.output(gc()))
    }
    else{
      plt <- plotAvgProf(tagMatrices, xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
        scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))) +
        scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)])))
      invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_peaks.png', sep=''), plot=plt, dpi=320)))
      invisible(capture.output(gc()))
    }
    plt <- plotAvgProf(tagMatrices[['Shared']], xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
      scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[length(conditions_colour_code)]) +
      scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[length(conditions_colour_code)])
    invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_shared-peaks.png', sep=''), plot=plt, dpi=320)))
    invisible(capture.output(gc()))
    rm(tagMatrices)
    rm(plt)
  },error = function(e)
  {
    message(e)
  }
)
invisible(capture.output(gc()))

# ========= Get Annotations =========
peakAnnoList <- list()
for (p in names(peaks)){
  cat("\nObtaining annotations for ", p, "\n")
  
  if (p == "Shared"){
    colour <- "#56B1F7"
    heat_colour <- "BuPu"
  }else{
    colour <- conditions_colour_code[[p]]
    heat_colour <- conditions_colour_code[[p]]
  }
  
  cat("\nAnnotating", p, '\n')
  anno <- annotatePeak(peaks[[p]], 
                       TxDb=anno_ref$txdb,
                       annoDb=anno_ref$annoDb,
                       level=opt$annotation_level,
                       tssRegion=c(-3000, 3000)
  )
  
  peakAnnoList[[p]] <- anno
  
  # Mapper for EntrezID to gene SYMBOL
  df <- as.data.frame(anno@anno)
  mapper <- df[c('geneId', 'SYMBOL')]
  mapper <- mapper[!duplicated(mapper), ]
  
  plt <- upsetplot(anno, vennpie=TRUE) + ggtitle(p)
  invisible(capture.output(ggsave(paste(result_dirs[[p]], p, '_annotated_peaks_upsetplot.png', sep=''), plot=plt, dpi=320, bg='white')))
  
  tryCatch(
    {
      png(paste(result_dirs[[p]], p, '_peaks_annotation_pie.png', sep=''), width=1680, height=1200)
      plt <- plotAnnoPie(anno, main=paste(p, '\n\n', length(anno@anno), ' Sites', sep=''), line=-10, cex.main=3.25, cex=3)
      invisible(capture.output( dev.off() ))
      #invisible(capture.output( ggsave(filename=paste(result_dirs[[p]], p, '_peaks_annotation_pie.png', sep=''), plot=grid.arrange(plt)) ))
    },error = function(e)
    {
      message(e)
    }
  )
  
  # Write annotation to file
  write.table(anno, file=paste(result_dirs[[p]], p, '_annotated.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
  
  tryCatch(
    {
      genes <- list()
      genes[[p]] <- anno@anno$geneId
      names(genes) = sub("_", "\n", names(genes))
      cat("\nGetting", p, 'KEGG\n')
      
      # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
      compKEGG <- compareCluster(geneCluster=genes,
                                 fun="enrichKEGG",
                                 pvalueCutoff=opt$fdr,
                                 pAdjustMethod="BH",
                                 minGSSize = opt$minGSSize,
                                 maxGSSize = opt$maxGSSize,
                                 organism=anno_ref$keggOrg
                                ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
      
      if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
        # Map EntrezIDs to gene SYMBOL
        compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
        myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
        for (i in 1:length(myEntrez)){
          compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
        }
        #plt <- dotplot(compKEGG, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
        plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=colour, n=15)
        invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], p, '_annotated_KEGG.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
        
        # Write annotations to csv
        write.table(as.data.frame(compKEGG), file=paste(result_dirs[[p]], p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
      } else{
        cat("\nNo annotation results\n")
      }
    },error = function(e)
    {
      message(e)
    }
  )
  invisible(capture.output(gc()))
  
  for (ont in c('ALL', 'CC', 'MF', 'BP')){
    tryCatch(
      {
        genes <- list()
        genes[[p]] <- anno@anno$SYMBOL
        names(genes) = sub("_", "\n", names(genes))
        cat("\nGetting ", p, ' GO ', ont, '\n')
        
        compGO <- compareCluster(geneCluster=genes,
                                 keyType='SYMBOL',
                                 OrgDb=anno_ref$annoDb,
                                 fun="enrichGO",
                                 ont=ont,
                                 pvalueCutoff=opt$fdr,
                                 pAdjustMethod="BH",
                                 minGSSize = opt$minGSSize,
                                 maxGSSize = opt$maxGSSize,
                                 readable=TRUE
        ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
        if ((!is.null(compGO)) & (dim(compGO@compareClusterResult)[1] > 0)){
          compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
          #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
          plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=colour, n=15)
          invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], p, '_annotated_GO-', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
          
          # Write annotations to csv
          write.table(as.data.frame(compGO), file=paste(result_dirs[[p]], p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
        } else{
        cat("\nNo annotation results\n")
        }
      },error = function(e)
      {
        message(e)
      }
    )
    invisible(capture.output(gc()))
    
  }

  genes_entrez <- list()
  genes_entrez[[p]] <- anno@anno$geneId
  names(genes_entrez) = sub("_", "\n", names(genes_entrez))
  for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
    cat("\nGetting ", p, ' DAVID ', annotation_type, ' annotations...\n')  
    # DAVID Annotation
    tryCatch(
      {
        compDAVID <- enrichDAVID(
                  unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))]),
                  idType = "ENTREZ_GENE_ID",
                  minGSSize = opt$minGSSize,
                  maxGSSize = opt$maxGSSize,
                  annotation = annotation_type,
                  pvalueCutoff = opt$fdr,
                  pAdjustMethod = "BH",
                  #species = NA,
                  david.user=opt$david_user
                )
        
        if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
          # Map EntrezIDs to gene SYMBOL
          compDAVID@result$SYMBOL <- compDAVID@result$geneID
          myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
          for (i in 1:length(myEntrez)){
            compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
          }
          plt <- make_dotplot(compDAVID@result, title=paste('DAVID - ', p, sep=""), ylabel=paste(annotation_type,"Category", sep=' '), colour=colour, n=15)
          invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], 'DAVID_annotation_', annotation_type, '_', p, '_dotplot.png', sep=''), plot=plt, dpi=320)))
          remove(plt)
          # Write annotations to csv
          write.table(as.data.frame(compDAVID), file=paste(result_dirs[[p]], 'DAVID_annotation_', annotation_type, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)

        } else{
        cat("\nNo annotation results\n")
        }
      },error = function(e)
      {
        message(e)
      }
    )
  }
  
}
          
tryCatch(
  {
    plt <- plotAnnoBar(peakAnnoList)
    invisible(capture.output(ggsave(filename=paste(result_dir, 'peaks_annotation_distribution_bar.png', sep=''), plot=plt, dpi=320)))
  },error = function(e)
  {
    message(e)
  }
)
tryCatch(
  {
    plt <- plotDistToTSS(peakAnnoList)
    invisible(capture.output(ggsave(filename=paste(result_dir, 'peaks_annotation_TSS_distribution.png', sep=''), plot=plt, dpi=320)))
  },error = function(e)
  {
    message(e)
  }
)
invisible(capture.output(gc()))

          
# ========= END OF OCCUPANCY ANALYSIS =========

# Free up memory
rm(dbObj.final)
rm(dbObj.total)
rm(plt)
rm(anno)
rm(peakAnnoList)
rm(compKEGG)
rm(compGO)
invisible(capture.output(gc()))


if (opt$occupancy_only){
  cat("Skipping affinity analysis.\nDone!\n")
  q()
}
          
# ========= START OF AFFINITY ANALYSIS =========
# The consensus peaks determined above are used here
# to focus on relevant peak regions only. This way, read
# counts for consensus peaks are analyzed for significance 
# in DE.
cat("\n=====  Affinity Analysis =====\n")
result_dir <- paste(opt$result_dir, 'Affinity_Analysis/', sep='')
# Make new directory for analysis
if (!file.exists(result_dir)) {
  dir.create(result_dir)
}

supplementary_dir <- paste(result_dir, "Supplementary/", sep='')
if (!file.exists(supplementary_dir)) {
  dir.create(supplementary_dir)
}

if (!file.exists(paste(result_dir, 'DESeq2/', sep=''))) {
  dir.create(paste(result_dir, 'DESeq2/', sep=''))
}
if (!file.exists(paste(result_dir, 'edgeR/', sep=''))) {
  dir.create(paste(result_dir, 'edgeR/', sep=''))
}


change_dirs <- function(res_dir, method, subfolder){
  if ((method == DBA_DESEQ2) | (method == 'DESeq2')){
    method_dir <- paste(res_dir, 'DESeq2/', sep='')
    if (!file.exists(paste(method_dir, subfolder, '/', sep=''))) {
      dir.create(paste(method_dir, subfolder, '/', sep=''))
    }
  }else if ((method == DBA_EDGER) | (method == 'edgeR')){
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

output_prefix <- change_dirs(result_dir, '', '')

# Count fragments for peaks from bam files
#dbObj.counted <- dba.count(dbObj.caller_consensus, bUseSummarizeOverlaps=TRUE, 
#                           peaks=consensus_peaks, 
#                           minOverlap=1, 
#                           score=DBA_SCORE_NORMALIZED,
#                           fragmentSize=dbObj.caller_consensus$config$fragmentSize,
#summits=200, filter=1, bRemoveDuplicates=FALSE, bScaleControl=TRUE,
#bSubControl=is.null(dbObj.noblacklist$greylist),
#mapQCth=dbObj.noblacklist$config$mapQCth, 
#filterFun=max, minCount=0,
#bLog=FALSE,
#readFormat=DBA_READS_DEFAULT, 
#bParallel=dbObj.noblacklist$config$RunParallel
#)
# ========= Count Reads from bams for All Relevant Consensus Peaks =========
#(recommend to use ‘summits=100‘ for ATAC-seq).
if (length(unique(dbObj$samples$Factor)) > 1){
  for (n in names(dbObj.caller_consensus$config)){
    if (n != 'fragmentSize'){
      dbObj.caller_consensus$config[[n]] <- unique(dbObj.caller_consensus$config[[n]])
    }
  }
  dbObj.counted <- dba.count(dbObj.caller_consensus, bUseSummarizeOverlaps=TRUE, 
                             peaks=consensus_peaks, 
                             minOverlap=1, 
                             score=DBA_SCORE_NORMALIZED,
                             bParallel=TRUE,
                             #fragmentSize=dbObj.caller_consensus$config$fragmentSize,
                             #summits=200, filter=1, bRemoveDuplicates=FALSE, bScaleControl=TRUE,
                             #bSubControl=is.null(dbObj.noblacklist$greylist),
                             #mapQCth=dbObj.noblacklist$config$mapQCth, 
                             #filterFun=max, minCount=0,
                             #bLog=FALSE,
                             #readFormat=DBA_READS_DEFAULT, 
                             #bParallel=dbObj.noblacklist$config$RunParallel
  )
}else{
  for (n in names(dbObj.noblacklist$config)){
    if (n != 'fragmentSize'){
      dbObj.noblacklist$config[[n]] <- unique(dbObj.noblacklist$config[[n]])
    }
  }
  dbObj.counted <- dba.count(dbObj.noblacklist, bUseSummarizeOverlaps=TRUE, 
                             peaks=consensus_peaks, 
                             minOverlap=1, 
                             score=DBA_SCORE_NORMALIZED,
                             bParallel=TRUE,
                             #fragmentSize=dbObj.noblacklist$config$fragmentSize,
                             #summits=200, filter=1, bRemoveDuplicates=FALSE, bScaleControl=TRUE,
                             #bSubControl=is.null(dbObj.noblacklist$greylist),
                             #mapQCth=dbObj.noblacklist$config$mapQCth, 
                             #filterFun=max, minCount=0,
                             #bLog=FALSE,
                             #readFormat=DBA_READS_DEFAULT, 
                             #bParallel=dbObj.noblacklist$config$RunParallel
  )
}
dbObj.counted

# ========= Normalize Counts =========
tryCatch(
  {
    # Normalize ("safest general method")
    dbObj.norm <- dba.normalize(dbObj.counted, method=DBA_ALL_METHODS, 
                                normalize=DBA_NORM_NATIVE,
                                background=TRUE, library=DBA_LIBSIZE_DEFAULT,
                                spikein=FALSE, offsets=FALSE,
                                libFun=mean, bRetrieve=FALSE)
  },error = function(e)
  { 
    message(e, "\nUsing approximate normalization...")
    # Approximate normalization to above without extra reading of bam files
    dbObj.norm <- dba.normalize(dbObj.counted, method=DBA_ALL_METHODS, 
                                normalize=DBA_NORM_LIB,
                                background=FALSE, library=DBA_LIBSIZE_FULL,
                                spikein=FALSE, offsets=FALSE,
                                libFun=mean, bRetrieve=FALSE)
  }
)
cat("After normalizing:\n")
dbObj.norm

# Plots
png(paste(supplementary_dir, 'consensus_peaks_counted_normalized_heatmap.png', sep=''))
dba.plotHeatmap(dbObj.norm)
invisible(capture.output( dev.off() ))

png(paste(supplementary_dir, 'consensus_peaks_counted_normalized_pca.png', sep=''))
dba.plotPCA(dbObj.norm, attributes=DBA_CONDITION, label=DBA_ID, vColors=(colours))
invisible(capture.output( dev.off() ))
invisible(capture.output(gc()))

# ========= Define contrasts between sample conditions (assumes 2 conditions) =========
# Define design and contrasts explicitly
if (length(unique(dbObj.norm$samples$Replicate)) > 1) {
  dbObj.contrast <- dba.contrast(dbObj.norm, design="~Condition", 
                                 minMembers=2, 
                                 contrast=c("Condition", unique(dbObj.norm$samples$Condition)))
}else{
  # Pre-DiffBind version 3 (or if no use of multiple replicates...will not be able to perform analysis in this case)
  dbObj.contrast <- dba.contrast(dbObj.norm, design=FALSE, categories=DBA_CONDITION,
                                 group1=dbObj.norm$masks[[unique(dbObj.norm$samples$Condition)[1]]], 
                                 name1=unique(dbObj.norm$samples$Condition)[1],
                                 name2=unique(dbObj.norm$samples$Condition)[2],
                                 minMembers=2,
                                 reorderMeta=list(Condition=unique(dbObj.norm$samples$Condition)[1], 
                                                  Condition=unique(dbObj.norm$samples$Condition)[2])
  )
}
cat("After defining contrast:\n")
dbObj.contrast

# ========= Analyze Affinities =========
dbObj.analyzed <- dba.analyze(dbObj.contrast, method=DBA_ALL_METHODS, bParallel=TRUE)
cat("After analyzing:\n")
dbObj.analyzed


# Plots
tryCatch(
  {
    png(paste(result_dir, 'analyzed_venn.png', sep=''))
    plt <- dba.plotVenn(main="DE Binding Sites Identified by Method",
                        dbObj.analyzed, 
                        contrast=1, 
                        method=DBA_ALL_METHODS
    )
    invisible(capture.output( dev.off() ))
    #invisible(capture.output( ggsave(filename=paste(result_dir, 'analyzed_venn.png', sep=''), plot=grid.arrange(plt)) ))
  }, error=function(e){
    message("No figure\n", e)
    if (file.exists(paste(result_dir, 'analyzed_venn.png', sep=''))){
      invisible(file.remove(paste(result_dir, 'analyzed_venn.png', sep='')))
    }
    
  }
)
invisible(capture.output(gc()))

for (m in c(DBA_DESEQ2, DBA_EDGER)){
  
  output_prefix <- change_dirs(result_dir, m, '')
  
  tryCatch(
    {
      tryCatch(
        {
          v <- dba.plotVenn(dbObj.analyzed, method=m, contrast=1, 
                            bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE,
                            main=paste("Binding Sites - ", m, sep=''))
          png(paste(output_prefix, 'analyzed_venn_', m, '.png', sep=''))
          grid.newpage()
          g = draw.pairwise.venn(area1=length(v$onlyA)+length(v$inAll), 
                                 area2=length(v$onlyB)+length(v$inAll),
                                 cross.area=length(v$inAll),
                                 category=names(unique_peaks),
                                 fill=unname(unlist(conditions_colour_code)),
                                 col=NA, cex=1.5, cat.cex=1.5,
                                 cat.dist = c(0,0))
          grid.arrange(gTree(children=g), top=textGrob("Binding Sites", gp=gpar(cex=1.5, fontface='bold')),
                       bottom=textGrob(m, gp=gpar(cex=1.5)))
          
          invisible(capture.output(dev.off()))
        }, error=function(e){
          message('Peaksets do not meet specified criteria for venn')
          
          tryCatch(
            {
              png(paste(output_prefix, 'analyzed_venn_', m, '.png', sep=''))
              v <- dba.plotVenn(dbObj.analyzed, method=m, contrast=1, 
                                bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=TRUE,
                                main=paste("Binding Sites - ", m, sep=''))
              invisible(capture.output(dev.off()))
            }, error=function(e){
              message('Peaksets do not meet specified criteria for venn')
            }
          )
        }
      )
      
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      if (file.exists(paste(output_prefix, 'analyzed_venn_', m, '.png', sep=''))){
        invisible(file.remove(paste(output_prefix, 'analyzed_venn_', m, '.png', sep='')))
      }
      
    }
  )
  invisible(capture.output(gc()))
  
  tryCatch(
    {
      png(paste(output_prefix, 'analyzed_volcano_', m, '.png', sep=''))
      par(mfrow=c(1,2))
      dba.plotVolcano(dbObj.analyzed, method=m)
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      if (file.exists(paste(output_prefix, 'analyzed_volcano_', m, '.png', sep=''))){
        invisible(file.remove(paste(output_prefix, 'analyzed_volcano_', m, '.png', sep='')))
      }
      
    }
  )
  invisible(capture.output(gc()))
  
  tryCatch(
    {
      png(paste(output_prefix, 'analyzed_ma_', m, '.png', sep=''))
      dba.plotMA(dbObj.analyzed, method=m)
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      if (file.exists(paste(output_prefix, 'analyzed_ma_', m, '.png', sep=''))){
        invisible(file.remove(paste(output_prefix, 'analyzed_ma_', m, '.png', sep='')))
      }
      
    }
  )
  invisible(capture.output(gc()))
  
  tryCatch(
    {
      png(paste(output_prefix, 'analyzed_box_', m, '.png', sep=''))
      dba.plotBox(dbObj.analyzed, method=m)
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      if (file.exists(paste(output_prefix, 'analyzed_box_', m, '.png', sep=''))){
        invisible(file.remove(paste(output_prefix, 'analyzed_box_', m, '.png', sep='')))
      }
      
    }
  )
  invisible(capture.output(gc()))
  
}

output_prefix <- change_dirs(result_dir, '', '')

tryCatch(
  {
    profile_colors <- list()
    for (i in 1:length(unique(dbObj.analyzed$samples$Condition))) {
      profile_colors[[unique(dbObj.analyzed$samples$Condition)[i]]] <- c('white', colours[[i]])
    }
    #profile_colors <- rev(profile_colors)
    profiles <- dba.plotProfile(dbObj.analyzed, merge=c(DBA_REPLICATE), normalize=TRUE)
    png(paste(output_prefix, 'analyzed_profile.png', sep=''))
    #grid.newpage()
    g <- dba.plotProfile(profiles, matrices_color=profile_colors, all_color_scales_equal=FALSE, 
                         decreasing=FALSE, 
                         #group_anno_color=unlist(unname(conditions_colour_code))
    )
    #grid.arrange(gTree(children=List(g)), top=textGrob("Heatmap-Profiles", gp=gpar(cex=1.5, fontface='bold')),
    #             bottom=textGrob(m, gp=gpar(cex=1.5)))
    invisible(capture.output(dev.off()))
  }, error=function(e){
    message("No figure\n", e)
    invisible(capture.output(dev.off()))
    if (file.exists(paste(output_prefix, 'analyzed_profile.png', sep=''))){
      invisible(file.remove(paste(output_prefix, 'analyzed_profile.png', sep='')))
    }
    
  }
)
invisible(capture.output(gc()))


# ========= Generate DE Analysis Reports =========
reports <- list()
reports[["DESeq2"]] <- dba.report(dbObj.analyzed, method=DBA_DESEQ2, contrast=1, th=1)
reports[["edgeR"]] <- dba.report(dbObj.analyzed, method=DBA_EDGER, contrast=1, th=1)

for (report in names(reports)){
  
  output_prefix <- change_dirs(result_dir, report, '')
  
  # Write complete report to file
  res <- as.data.frame(reports[[report]])
  names(res)[names(res) == 'Fold'] <- 'log2FoldChange'
  names(res)[names(res) == 'FDR'] <- 'p.adjust'
  res <- res[order(res$log2FoldChange, decreasing=TRUE), ]
  res <- as.data.frame(annotatePeak(GRanges(res), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                    level=opt$annotation_level,
                                    tssRegion=c(-3000, 3000))@anno)
  
  write.table(res, file=paste(output_prefix, 'analyzed_report_', report, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
  
  if (is.null(dba.report(dbObj.analyzed, method=report, contrast=1, th=opt$fdr))){
    cat("\nNo DEGs identified by", report, "at a significance threshold of", opt$fdr, "skipping further analysis...\n")
    next
  }
  
  # Report columns are seqnames, start, end, width, strand, Conc, Conc_Group1, Conc_Group2, Fold, p.value, FDR
  # Create bed files for each keeping only significant peaks (p<0.05)
  # Comparing those whereby Fold > 0 vs Fold < 0, indicating enrichment gain vs loss of group 1 over group 2
  
  out <- as.data.frame(reports[[report]])
  #out <- out[as.data.frame(findOverlaps(both, reports[[report]]))[["subjectHits"]], ]
  
  cat("\n", report, "Report:\n")
  print(head(out))
  print(tail(out))
  
  gained <- out %>% 
    dplyr::filter(FDR < opt$fdr & Fold > opt$lfc) %>% 
    dplyr::select(seqnames, start, end)
  lost <- out %>% 
    dplyr::filter(FDR < opt$fdr & Fold < (0 - opt$lfc)) %>% 
    dplyr::select(seqnames, start, end)
  
  cat('\n', dim(out)[[1]], 'peaks in', report, 'report:\n')
  cat('\n\t', dim(gained)[[1]] + dim(lost)[[1]], 'statistically significant peaks', '( FDR <', dbObj.analyzed$config$th, ')\n')
  cat('\n\t', dim(gained)[[1]], 'DE peaks in', dbObj.contrast$contrasts[[1]]$name1, '(log2Fold-change > ', opt$lfc, '\n')
  cat('\n\t', dim(lost)[[1]], 'DE peaks in', dbObj.contrast$contrasts[[1]]$name2, '(log2Fold-change < -', opt$lfc, '\n')
  
  
  # Write DE result to bed files
  if (dim(gained)[1] > 0){
    output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name1)
    write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  }
  if (dim(lost)[1] > 0){
    output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
    write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  }
  
  # Plot profile heatmaps for all significant (FDR < opt$fdr) sites for each method
  output_prefix <- change_dirs(result_dir, report, '')
  tryCatch(
    {
      profile_colors <- list()
      for (i in 1:length(unique(dbObj.analyzed$samples$Condition))) {
        profile_colors[[unique(dbObj.analyzed$samples$Condition)[i]]] <- c('white', colours[[i]])
      }
      #profile_colors <- rev(profile_colors)
      # Plots all significant sites among both conditions (will include signals)
      profiles_significant <- dba.plotProfile(dbObj.analyzed, merge=c(DBA_REPLICATE), normalize=TRUE, 
                                              sites=GRanges(out %>% dplyr::filter(FDR < opt$fdr))
      )
      png(paste(output_prefix, 'analyzed_profiles_', report,'.png', sep=''))
      dba.plotProfile(profiles_significant, matrices_color=profile_colors, all_color_scales_equal=FALSE, 
                      decreasing=FALSE, 
                      #group_anno_color=unlist(unname(conditions_colour_code))
      )
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      if (file.exists(paste(output_prefix, 'analyzed_profiles_', report,'.png', sep=''))){
        invisible(file.remove(paste(output_prefix, 'analyzed_profiles_', report,'.png', sep='')))
      }
      
    }
  )
  invisible(capture.output(gc()))
  
  gained <- out %>% 
    dplyr::filter(FDR < opt$fdr & Fold > opt$lfc)
  lost <- out %>% 
    dplyr::filter(FDR < opt$fdr & Fold < (0 - opt$lfc))
  
  # Write complete DE result to files
  if (dim(gained)[1] > 0){
    output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name1)
    gained <- as.data.frame(annotatePeak(GRanges(gained), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                         level=opt$annotation_level,
                                         tssRegion=c(-3000, 3000))@anno)
    gained <- gained[order(gained$Fold, decreasing=TRUE), ]
    write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
  }
  if (dim(lost)[1] > 0){
    output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
    lost <- as.data.frame(annotatePeak(GRanges(lost), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                       level=opt$annotation_level,
                                       tssRegion=c(-3000, 3000))@anno)
    lost <- lost[order(lost$Fold, decreasing=FALSE), ]
    write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
  }
  
  # Use absolute fold-change for plotting magnitude
  gained$Fold <- abs(gained$Fold)
  lost$Fold <- abs(lost$Fold)
  
  gpeaks <- GenomicRanges::GRangesList(Gained=gained, Lost=lost)
  names(gpeaks) <- c(dbObj.contrast$contrasts[[1]]$name1, dbObj.contrast$contrasts[[1]]$name2)
  
  # Plot peaks gained/lost over genome between conditions
  output_prefix <- change_dirs(result_dir, report, '')
  tryCatch(
    {
      plt <- covplot(gpeaks, title=paste("Peaks over Genome", report, sep=' - '), weightCol='Fold') + 
        scale_color_manual(values=rev(c(colours[1:length(unique_peaks)]))) + 
        scale_fill_manual(values=rev(c(colours[1:length(unique_peaks)])))
      invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_significant_merged_peaks.png', sep=''), plot=plt, dpi=320)))
      plt <- plt + facet_grid(chr ~ .id)
      invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_significant_peaks.png', sep=''), plot=plt, dpi=320)))
    }, error=function(e){
      message("No figure\n", e)
    }
  )
  invisible(capture.output(gc()))
  
  # ========= Get Annotations =========
  peakAnnoList <- list()
  for (p in names(gpeaks)){
    cat("\nObtaining annotations for ", p, "\n")
    output_prefix <- change_dirs(result_dir, report, p)
    
    if (length(gpeaks[[p]]) == 0){
      cat("\nNo", p, "peaks to annotate...\n")
      next
    }
    
    if (p == "Shared"){
      colour <- "#56B1F7"
      heat_colour <- "BuPu"
    }else{
      colour <- conditions_colour_code[[p]]
    }
    
    if (colour == "#00BFC4"){
      heat_colour <- "Blues"
    }else if (colour == "#F8766D"){
      heat_colour <- "Reds"
    }else if (colour == "#7CAE00"){
      heat_colour <- "Greens"
    }else if (colour == "#C77CFF"){
      heat_colour <- "BuPu"
    }else if (colour == "#e69e02"){
      heat_colour <- "YlOrBr"
    }else{
      heat_colour <- "Greys"
    }
    
    if (length(gpeaks[[p]]) > 0){
      cat("\nAnnotating", p, "\n")
      anno <- annotatePeak(gpeaks[[p]], 
                           TxDb=anno_ref$txdb,
                           annoDb=anno_ref$annoDb,
                           level=opt$annotation_level,
                           tssRegion=c(-3000, 3000))
      peakAnnoList[[p]] <- anno
      
      tryCatch(
        {
          png(paste(output_prefix,'_', report, 'peaks_annotation_pie', p, '.png', sep=''), width=1680, height=1200)
          plt <- plotAnnoPie(anno, main=paste(p, '\n\n', length(anno@anno), ' Sites', sep=''), line=-10, cex.main=3.25, cex=3)
          invisible(capture.output( dev.off() ))
          #invisible(capture.output( ggsave(filename=paste(output_prefix, report, 'peaks_annotation_pie', p, '.png', sep=''), plot=grid.arrange(plt)) ))
        },error = function(e)
        {
          message(e)
        }
      )
      
      plt <- upsetplot(anno, vennpie=TRUE) + ggtitle(p)
      invisible(capture.output(ggsave(paste(output_prefix, report, '_', p, '_annotated_peaks_upsetplot.png', sep=''), plot=plt, dpi=320, bg='white')))
      
      # Write annotation to file
      write.table(anno, file=paste(output_prefix, 'annotated_', report, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
      
      tryCatch(
        {
          genes <- list()
          genes[[p]] <- anno@anno$geneId
          names(genes) = sub("_", "\n", names(genes))
          cat("\nGetting ", report, ' ', p, ' KEGG\n')
          
          # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
          compKEGG <- compareCluster(geneCluster=genes,
                                     fun="enrichKEGG",
                                     pvalueCutoff=opt$fdr,
                                     pAdjustMethod="BH",
                                     minGSSize = opt$minGSSize,
                                     maxGSSize = opt$maxGSSize,
                                     organism=anno_ref$keggOrg
                                    ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          
          if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
            # Map EntrezIDs to gene SYMBOL
            compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
            myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
            for (i in 1:length(myEntrez)){
              compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
            }
            
            plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=colour, n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_KEGG.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
            
            # Write annotations to csv
            write.table(as.data.frame(compKEGG), file=paste(output_prefix, report, '_', p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
            
            plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, ylabel="KEGG Category")
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, ylabel="KEGG Category")
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            
          } else{
            cat("\nNo annotation results\n")
          }
        },error = function(e)
        {
          message(e)
        }
      )
      invisible(capture.output(gc()))
      
      for (ont in c('ALL', 'CC', 'MF', 'BP')){
        cat("\nGetting ", report, ' ',  p, ' GO - ', ont, '\n')
        
        tryCatch(
          {
            compGO <- compareCluster(geneCluster=genes,
                                     OrgDb=anno_ref$annoDb,
                                     fun="enrichGO",
                                     ont=ont,
                                     pvalueCutoff=opt$fdr,
                                     pAdjustMethod="BH",
                                     minGSSize = opt$minGSSize,
                                     maxGSSize = opt$maxGSSize,
                                     readable=TRUE
            ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
            
            if (!is.null(compGO)){
              compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
              #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
              plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=colour, n=15)
              invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
              
              # Write annotations to csv
              write.table(as.data.frame(compGO), file=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_pheatmapplot(compGO@compareClusterResult, res, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
              invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap_by_gene.png', sep=''), plot=plt, dpi=320)))
              remove(plt)
              plt <- make_pheatmapplot(compGO@compareClusterResult, res, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
              invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
              remove(plt)
            } else{
              cat("\nNo annotation results\n")
            }
          },error = function(e)
          {
            message(e)
          }
        )
        invisible(capture.output(gc()))
        
      }
        genes_entrez <- list()
        genes_entrez[[p]] <- anno@anno$geneId
        names(genes_entrez) = sub("_", "\n", names(genes_entrez))
        for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
          cat("\nGetting ", report, ' ', p, ' DAVID ', annotation_type, ' annotations...\n')  
          # DAVID Annotation
          tryCatch(
            {
              compDAVID <- enrichDAVID(
                        unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))]),
                        idType = "ENTREZ_GENE_ID",
                        minGSSize = opt$minGSSize,
                        maxGSSize = opt$maxGSSize,
                        annotation = annotation_type,
                        pvalueCutoff = opt$fdr,
                        pAdjustMethod = "BH",
                        #species = NA,
                        david.user=opt$david_user
                      )
              
              if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
                # Map EntrezIDs to gene SYMBOL
                compDAVID@result$SYMBOL <- compDAVID@result$geneID
                myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
                for (i in 1:length(myEntrez)){
                  compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                }
                
                plt <- make_dotplot(compDAVID@result, title=paste('DAVID - ', p, sep=""), ylabel=paste(annotation_type,"Category", sep=' '), colour=colour, n=15)
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '_dotplot.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                # Write annotations to csv
                write.table(as.data.frame(compDAVID), file=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_pheatmapplot(compDAVID@result, res, anno_type="DAVID", assembly=opt$assembly, title=paste('DAVID - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE)
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                plt <- make_pheatmapplot(compDAVID@result, res, anno_type="DAVID", assembly=opt$assembly, title=paste('DAVID - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE)
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
              } else{
                cat("\nNo annotation results\n")
              }
            },error = function(e)
            {
              message(e)
            }
          )
        }
    }
    else{
      cat("\nNo peaks to annotate for", p, '\n')
    }
  }
  
  output_prefix <- change_dirs(result_dir, report, '')
  tryCatch(
    {
      plt <- plotAnnoBar(peakAnnoList)
      invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))
    },error = function(e)
    {
      message(e)
    }
  )
  tryCatch(
    {    
      plt <- plotDistToTSS(peakAnnoList)
      invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_TSS_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))
    },error = function(e)
    {
      message(e)
    }
  )
  invisible(capture.output(gc()))
}

# Quit if one method produced no results to avoid repeating functions...
if (is.null(dba.report(dbObj.analyzed, method=DBA_DESEQ2, contrast=1, th=opt$fdr)) | is.null(dba.report(dbObj.analyzed, method=DBA_EDGER, contrast=1, th=opt$fdr))){
  cat("\nNo DEGs identified by either DESeq2 or edgeR at a significance threshold of", opt$fdr, "skipping further analysis...\n")
  cat("\nFINISHED!\n")
  q()
}


# For considering sites identified by DESeq2 AND edgeR
#reports[["DESeq2_and_edgeR"]] <- c(reports[["DESeq2"]], reports[["edgeR"]])
#names(reports[["DESeq2_and_edgeR"]]) <- 1:length(names(reports[["DESeq2_and_edgeR"]]))
df_both <- merge(as.data.frame(reports[["DESeq2"]]), as.data.frame(reports[["edgeR"]]), by=c("seqnames", "start", "end"), suffixes=c(".DESeq2", ".edgeR"))
reports[["DESeq2_and_edgeR"]] <- GRanges(df_both)
report <- "DESeq2_and_edgeR"

output_prefix <- change_dirs(result_dir, report, '')

# Write complete report to file
res <- as.data.frame(reports[[report]])
res['log2FoldChange'] <- (res$Fold.DESeq2 + res$Fold.edgeR)/2
res['p.adjust'] <- (res$FDR.DESeq2 + res$FDR.edgeR)/2
res <- as.data.frame(annotatePeak(GRanges(res), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                  level=opt$annotation_level,
                                  tssRegion=c(-3000, 3000))@anno)
res <- res[order(res$log2FoldChange, decreasing=TRUE), ]

write.table(res, file=paste(output_prefix, 'analyzed_report_', report, '.tsv', sep=''), sep="\t", quote=F, row.names=F)

out <- as.data.frame(reports[[report]])
#out <- out[as.data.frame(findOverlaps(both, reports[[report]]))[["subjectHits"]], ]

cat("\n", report, "Report:\n")
print(head(out))
print(tail(out))


gained <- out %>% 
  dplyr::filter(FDR.DESeq2 < opt$fdr & Fold.DESeq2 > opt$lfc & FDR.edgeR < opt$fdr & Fold.edgeR > opt$lfc) %>% 
  dplyr::select(seqnames, start, end)
lost <- out %>% 
  dplyr::filter(FDR.DESeq2 < opt$fdr & Fold.DESeq2 < (0 - opt$lfc) & FDR.edgeR < opt$fdr & Fold.edgeR < (0 - opt$lfc)) %>% 
  dplyr::select(seqnames, start, end)

cat('\n', dim(out)[[1]], 'peaks in', report, 'report:\n')
cat('\n\t', dim(gained)[[1]] + dim(lost)[[1]], 'statistically significant peaks', '( FDR <', dbObj.analyzed$config$th, ')\n')
cat('\n\t', dim(gained)[[1]], 'DE peaks log2FoldChange >', opt$lfc, 'in', dbObj.contrast$contrasts[[1]]$name1, '\n')
cat('\n\t', dim(lost)[[1]], 'DE peaks log2FoldChange <', (0 - opt$lfc), 'in', dbObj.contrast$contrasts[[1]]$name2, '\n')


# Write DE result to bed files
if (dim(gained)[1] > 0){
  output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name1)
  write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
}
if (dim(lost)[1] > 0){
  output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
  write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
}

gained <- out %>% 
  dplyr::filter(FDR.DESeq2 < opt$fdr & Fold.DESeq2 > opt$lfc & FDR.edgeR < opt$fdr & Fold.edgeR > opt$lfc)
lost <- out %>% 
  dplyr::filter(FDR.DESeq2 < opt$fdr & Fold.DESeq2 < (0 - opt$lfc) & FDR.edgeR < opt$fdr & Fold.edgeR < (0 - opt$lfc))


if (dim(gained)[1] > 0){
  output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name1)
  gained <- as.data.frame(annotatePeak(GRanges(gained), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                       level=opt$annotation_level,
                                       tssRegion=c(-3000, 3000))@anno)
  gained <- gained[order(gained$Fold, decreasing=TRUE), ]
  write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
}
if (dim(lost)[1] > 0){
  output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
  lost <- as.data.frame(annotatePeak(GRanges(lost), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                     level=opt$annotation_level,
                                     tssRegion=c(-3000, 3000))@anno)
  lost <- lost[order(lost$Fold, decreasing=FALSE), ]
  write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
}

# Use absolute fold-change for plotting magnitude
#gained$Fold <- abs(gained$Fold)
#lost$Fold <- abs(lost$Fold)

gpeaks <- GenomicRanges::GRangesList(Gained=gained, Lost=lost)
names(gpeaks) <- c(dbObj.contrast$contrasts[[1]]$name1, dbObj.contrast$contrasts[[1]]$name2)

# Plot peaks gained/lost over genome between conditions
output_prefix <- change_dirs(result_dir, report, '')
tryCatch(
  {
    plt <- covplot(gpeaks, title=paste("Peaks over Genome", report, sep=' - '), ) + #weightCol='Fold') + 
      scale_color_manual(values=rev(c(colours[1:length(unique_peaks)]))) + 
      scale_fill_manual(values=rev(c(colours[1:length(unique_peaks)])))
    invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_significant_merged_peaks.png', sep=''), plot=plt, dpi=320)))
    plt <- plt + facet_grid(chr ~ .id)
    invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_significant_peaks.png', sep=''), plot=plt, dpi=320)))
  }, error=function(e){
    message("No figure\n", e)
    invisible(capture.output(dev.off()))
  }
)
invisible(capture.output(gc()))

# ========= Get Annotations =========
peakAnnoList <- list()
for (p in names(gpeaks)){
  cat("\nObtaining annotations for ", p, "\n")
  output_prefix <- change_dirs(result_dir, report, p)
  
  if (length(gpeaks[[p]]) == 0){
    cat("\nNo", p, "peaks to annotate...\n")
    next
  }
  
  if (p == "Shared"){
    colour <- "#56B1F7"
    heat_colour <- "BuPu"
  }else{
    colour <- conditions_colour_code[[p]]
    heat_colour <- conditions_colour_code[[p]]
  }
  
  if (colour == "#00BFC4"){
    heat_colour <- "Blues"
  }else if (colour == "#F8766D"){
    heat_colour <- "Reds"
  }else if (colour == "#7CAE00"){
    heat_colour <- "Greens"
  }else if (colour == "#C77CFF"){
    heat_colour <- "BuPu"
  }else if (colour == "#e69e02"){
    heat_colour <- "YlOrBr"
  }else{
    heat_colour <- "Greys"
  }
  
  if (length(gpeaks[[p]]) > 0){
    cat("\nAnnotating", p, "\n")
    anno <- annotatePeak(gpeaks[[p]], 
                         TxDb=anno_ref$txdb,
                         annoDb=anno_ref$annoDb,
                         level=opt$annotation_level,
                         tssRegion=c(-3000, 3000))
    peakAnnoList[[p]] <- anno
    
    tryCatch(
      {
        png(paste(output_prefix,'_', report, 'peaks_annotation_pie', p, '.png', sep=''), width=1680, height=1200)
        plt <- plotAnnoPie(anno, main=paste(p, '\n\n', length(anno@anno), ' Sites', sep=''), line=-10, cex.main=3.25, cex=3)
        invisible(capture.output( dev.off() ))
        #invisible(capture.output( ggsave(filename=paste(output_prefix, report, 'peaks_annotation_pie', p, '.png', sep=''), plot=grid.arrange(plt)) ))
      },error = function(e)
      {
        message(e)
      }
    )
    
    plt <- upsetplot(anno, vennpie=TRUE) + ggtitle(p)
    invisible(capture.output(ggsave(paste(output_prefix, report, '_', p, '_annotated_peaks_upsetplot.png', sep=''), plot=plt, dpi=320, bg='white')))
    
    # Write annotation to file
    write.table(anno, file=paste(output_prefix, 'annotated_', report, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
    
    tryCatch(
      {
        genes <- list()
        genes[[p]] <- anno@anno$geneId
        names(genes) = sub("_", "\n", names(genes))
        cat("\nGetting", report, p, 'KEGG\n')
        
        # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
        compKEGG <- compareCluster(geneCluster=genes,
                                   fun="enrichKEGG",
                                   pvalueCutoff=opt$fdr,
                                   pAdjustMethod="BH",
                                   minGSSize = opt$minGSSize,
                                   maxGSSize = opt$maxGSSize,
                                   organism=anno_ref$keggOrg
                                  ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
        
        if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
          # Map EntrezIDs to gene SYMBOL
          compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
          myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
          for (i in 1:length(myEntrez)){
            compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
          }
          
          plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=colour, n=15)
          invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_KEGG.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
          
          # Write annotations to csv
          write.table(as.data.frame(compKEGG), file=paste(output_prefix, report, '_', p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
          
          plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, ylabel="KEGG Category")
          invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
          remove(plt)
          plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, ylabel="KEGG Category")
          invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
          remove(plt)
        } else{
          cat("\nNo annotation results\n")
        }
      },error = function(e)
      {
        message(e)
      }
    )
    invisible(capture.output(gc()))
    
    for (ont in c('ALL', 'CC', 'MF', 'BP')){
      cat("\nGetting ", report, ' ',  p, ' GO - ', ont, '\n')
      tryCatch(
        {
          compGO <- compareCluster(geneCluster=genes,
                                   OrgDb=anno_ref$annoDb,
                                   fun="enrichGO",
                                   ont=ont,
                                   pvalueCutoff=opt$fdr,
                                   pAdjustMethod="BH",
                                   minGSSize = opt$minGSSize,
                                   maxGSSize = opt$maxGSSize,
                                   readable=TRUE
          ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          
          if ((!is.null(compGO)) & (dim(compGO@compareClusterResult)[1] > 0)){
            compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
            #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
            plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=colour, n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
            
            # Write annotations to csv
            write.table(as.data.frame(compGO), file=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
            
            plt <- make_pheatmapplot(compGO@compareClusterResult, res, assembly=opt$assembly, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap_by_gene.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            plt <- make_pheatmapplot(compGO@compareClusterResult, res, assembly=opt$assembly, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            
          } else{
            cat("\nNo annotation results\n")
          }
        },error = function(e)
        {
          message(e)
        }
      )
      invisible(capture.output(gc()))
      
    }
        
    genes_entrez <- list()
    genes_entrez[[p]] <- anno@anno$geneId
    names(genes_entrez) = sub("_", "\n", names(genes_entrez))
    for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
      cat("\nGetting ", report, ' ', p, ' DAVID ', annotation_type, ' annotations...\n')  
      # DAVID Annotation
      tryCatch(
        {
          compDAVID <- enrichDAVID(
                    unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))]),
                    idType = "ENTREZ_GENE_ID",
                    minGSSize = opt$minGSSize,
                    maxGSSize = opt$maxGSSize,
                    annotation = annotation_type,
                    pvalueCutoff = opt$fdr,
                    pAdjustMethod = "BH",
                    #species = NA,
                    david.user=opt$david_user
                  )
          
          if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
            # Map EntrezIDs to gene SYMBOL
            compDAVID@result$SYMBOL <- compDAVID@result$geneID
            myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
            for (i in 1:length(myEntrez)){
              compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
            }
            
            plt <- make_dotplot(compDAVID@result, title=paste('DAVID - ', p, sep=""), ylabel=paste(annotation_type,"Category", sep=' '), colour=colour, n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '_dotplot.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            # Write annotations to csv
            write.table(as.data.frame(compDAVID), file=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
            
            plt <- make_pheatmapplot(compDAVID@result, res, anno_type="DAVID", assembly=opt$assembly, title=paste('DAVID - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=TRUE)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
            plt <- make_pheatmapplot(compDAVID@result, res, anno_type="DAVID", assembly=opt$assembly, title=paste('DAVID - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=opt$lfc, dendro=TRUE, sort_genes=FALSE)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, 'DAVID_annotation_', annotation_type, '_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
            remove(plt)
          } else{
            cat("\nNo annotation results\n")
          }
        },error = function(e)
        {
          message(e)
        }
      )
    }
    
  }
  else{
    cat("No peaks to annotate for", p, '\n')
  }
}

output_prefix <- change_dirs(result_dir, report, '')
tryCatch(
  {
    plt <- plotAnnoBar(peakAnnoList)
    invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))
  },error = function(e)
  {
    message(e)
  }
)
tryCatch(
  {
    plt <- plotDistToTSS(peakAnnoList)
    invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_TSS_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))
  },error = function(e)
  {
    message(e)
  }
)
invisible(capture.output(gc()))

cat('\nFINISHED!\n')
