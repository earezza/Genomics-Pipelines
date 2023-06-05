# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
  library(DiffBind)
  library(tidyverse)
  #library(profileplyr)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  #library(ensembldb)
  library(EnsDb.Hsapiens.v86)
  library(EnsDb.Mmusculus.v79)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  #library(reshape2)
  library(ReactomePA)
  library(clusterProfiler)
  library(vulcan)
  library(UpSetR)
  library(dplyr)
  library(VennDiagram)
  library(eulerr)
  require(gridExtra)
  library(optparse)
}))

# ======= Get command-line optional arguments =======
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="DiffBind-formatted sample sheet", metavar="character"),
  make_option(c("-s", "--fragmentsizes"), type="character", default=NULL, help="File with bam fragment sizes generated from bamPEFragmentSize, otherwise determined via vulcan", metavar="character"),
  make_option(c("-o", "--organism"), type="character", default="mouse", help="Organism to annotate genes at peaks, human or mouse (default)", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="Peaks_Analysis/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-d", "--database"), type="character", default="ucsc", help="Database reference for peaks gene annotations, ucsc (default) or ensembl", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cat("Command-line options:\n")
for (i in which(names(opt) != "help")) {
  cat(names(opt)[i], '=', paste(opt)[i], "\n")
}

# ========= SETUP RUN AND VARIABLES =========
# Get samplesheet, set directory with samplesheet as workspace
samplesheet <- opt$file

# Get mean fragmentSize for samples from deepTools bamPEFragmentSize function (obtained prior to Rscript)
# Note: Assumes order of samples matches the samplesheet
if (is.null(opt$fragmentsizes)){
  fragment_size <- 125 # default
} else{
  fragment_sizes <- read.table(opt$fragmentsizes, header=TRUE, sep='\t')[,c( "X", "Frag..Len..Mean")]
  fragment_size <- round(fragment_sizes["Frag..Len..Mean"], 0)[,1]
  fragment_size <- round(mean(fragment_size))
}

# ========= Get database references for annotations =========
if (opt$organism == "mouse"){
  annoDb <- "org.Mm.eg.db"
  keggOrg <- "mmu"
  if (opt$database == "ucsc"){
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  }else if (opt$database == "ensemble"){
    txdb <- EnsDb.Mmusculus.v79
    seqlevelsStyle(txdb) <- "UCSC" # format ensembl genes using UCSC style
  }
  else{
    stop("Invalid choice of annotation database")
  }
} else if (opt$organism == "human"){
  annoDb <- "org.Hs.eg.db"
  keggOrg <- "hsa"
  if (opt$database == "ucsc"){
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }else if (opt$database == "ensemble"){
    txdb <- EnsDb.Hsapiens.v86
    seqlevelsStyle(txdb) <- "UCSC" # format ensembl genes using UCSC style
  }
  else{
    stop("Invalid choice of annotation database")
  }
} else{
  stop("Invalid choice of organism")
}

# Get promoter regions from database
promoters <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

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

# # ========= START OCCUPANCY ANALYSIS =========
# Here, peaks declared by peak caller(s) are used to identify
# differential expression between conditions. Overlapping peaks 
# between replicates/conditions are determined by the range of the peaks. 
# In this case, using occupancy alone for DE provides a less conservative
# analysis for DE by simply considering where peaks exist.
#
# Following this, affinity analysis below can provide a more conservative
# result for DE since the read counts are accounted for (peak shapes).


# ========= Load peaksets =========
setwd(dirname(samplesheet))
samplesheet <- basename(samplesheet)

# Set directories
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}
cat("Output files will be in", opt$result_dir, "\n")

result_dir <- paste(opt$result_dir, "Occupancy_Analysis/", sep='')
if (!file.exists(result_dir)) {
  dir.create(result_dir)
}

output_prefix <- gsub('.csv', '_', paste(result_dir, samplesheet, sep=""))
output_prefix <- gsub("diffbind_samplesheet_", "", output_prefix)

# Average overl all samples
#tryCatch(
#  {
#    png(paste(output_prefix, 'fragment_length.png', sep=""))
#    fragment_size <- average_fragment_length(unique(read.csv(samplesheet)$bamReads), plot=TRUE)
#    invisible(capture.output(dev.off()))
#  }, error=function(e){
#    message("Problem reading samplesheet.\n", e)
#    invisible(capture.output(dev.off()))
#    invisible(file.remove(paste(output_prefix, 'fragment_length.png', sep="")))
#  }
#)

# Average for each sample
fragment_size <- 1:length(read.csv(samplesheet)$SampleID)
for (b in unique(read.csv(samplesheet)$Condition)){
  for (r in unique(read.csv(samplesheet)$Replicate)){
    png(paste(output_prefix, 'fragment_length_', b, '-', r, '.png', sep=""))
    mean_fragment_size <- average_fragment_length(read.csv(samplesheet)$bamReads[[which(read.csv(samplesheet)$Condition == b)[1]]], plot=TRUE)
    for (i in which(read.csv(samplesheet)$Condition == b & read.csv(samplesheet)$Replicate == r)){
      fragment_size[i] <- mean_fragment_size
    }
    invisible(capture.output(dev.off())) 
  }
}
fragment_size <- 125 # default

dbObj <- dba(sampleSheet=samplesheet, minOverlap=1,
             config=data.frame(th=0.05,
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
colours <- c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF", "#e69e02")
conditions_colour_code <- list()
for (i in 1:length(unique(dbObj$samples$Condition))) {
  conditions_colour_code[[unique(dbObj$samples$Condition)[i]]] <- colours[i]
}

png(paste(output_prefix, 'raw_heatmap.png', sep=""))
dba.plotHeatmap(dbObj)
invisible(capture.output(dev.off()))

# Show overlap rates for each condition
cat("Peak overlaps in at least (1, 2, ...) replicates/callers for each condition:\n")
png(paste(output_prefix, "raw_overlap_rates.png", sep=""))
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

# ========= Remove Blacklisted Regions =========
# Remove blacklisted regions to ignore irrelevant peaks (blacklisted regions from ENCODE, genome selected is based on prediction from bam files)
dbObj.noblacklist <- dba.blacklist(dbObj, blacklist=TRUE, greylist=FALSE)
blacklisted_peaks <- dba.blacklist(dbObj.noblacklist, Retrieve=DBA_BLACKLISTED_PEAKS)
cat("After blacklist applied:\n")
dbObj.noblacklist

png(paste(output_prefix, 'raw_noblacklist_heatmap.png', sep=""))
dba.plotHeatmap(dbObj.noblacklist)
invisible(capture.output(dev.off()))

# Show overlap rates for each condition
cat("Peak overlaps in at least (1, 2, ...) replicates/callers for each condition:\n")
png(paste(output_prefix, "raw_noblacklist_overlap_rates.png", sep=""))
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

# ========= Get Consensus Peaks =========
# If using peaks from multiple peak callers (defined in Factor column of samplesheet)
if (length(unique(dbObj$samples$Factor)) > 1){
  # Get consensus overlaps in peaksets for each condition (peaks must be overlapping in 2/3rds of peak callers)
  dbObj.total <- dba.peakset(dbObj.noblacklist, consensus=c(DBA_CONDITION, DBA_REPLICATE), minOverlap=0.66)
  # resulting consensus between callers
  dbObj.caller_consensus <- dba(dbObj.total, mask=dbObj.total$masks$Consensus, minOverlap=1)
  if (length(unique(dbObj$samples$Replicate)) > 1){
    # Get consensus between replicates for each condition (peaks must be overlapping in at least 2 replicates)
    dbObj.final <- dba.peakset(dbObj.caller_consensus, consensus=c(DBA_CONDITION), minOverlap=2)
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
  # consensus between replicates
  dbObj.total <- dba.peakset(dbObj.noblacklist, consensus=c(DBA_CONDITION), minOverlap=2)
  maskname <- names(dbObj.total$masks)[grepl('Replicate.1-2', names(dbObj.total$masks))]
  if (length(maskname) == 1){
    dbObj.consensus <- dba(dbObj.total, mask=dbObj.total$masks[[maskname]], minOverlap=1)
  }else if (length(maskname) == 2){
    dbObj.consensus <- dba(dbObj.total, mask=(dbObj.total$masks[[maskname[1]]] | dbObj.total$masks[[maskname[2]]]), minOverlap=1)
  }
}
dbObj.consensus

# Re-sort colours if condition orders changed after consensus (occurs when one condition has only 1 replicate, consensus must be added manually...)
# dbObj.consensus <- dba(dbObj.final, mask=(dbObj.final$masks$`Replicate.1-2` | dbObj.final$masks$CONDITION_WITH_ONE_REPLICATE), minOverlap=1)
temp <- list()
for (i in 1:length(conditions_colour_code)){
  temp[names(conditions_colour_code[which(names(conditions_colour_code) == dba.show(dbObj.consensus)$Condition[i])])] <- conditions_colour_code[which(names(conditions_colour_code) == dba.show(dbObj.consensus)$Condition[i])]
}
conditions_colour_code <- temp

# Consensus peaks from all conditions (all relevant peaks)
consensus_peaks <- dba.peakset(dbObj.consensus, bRetrieve=TRUE)

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
cat("\n", length(shared_peaks[["Shared"]]), "shared peaks.\n")

# Plots
if (length(unique(dba.show(dbObj.consensus)$Condition)) == 2){
  png(paste(output_prefix, 'raw_consensus_peaks.png', sep=""))
  grid.newpage()
  g = draw.pairwise.venn(area1=length(unique_peaks[[names(unique_peaks)[1]]])+length(shared_peaks[[names(shared_peaks)[1]]]), 
                         area2=length(unique_peaks[[names(unique_peaks)[2]]])+length(shared_peaks[[names(shared_peaks)[1]]]),
                         cross.area=length(shared_peaks[[names(shared_peaks)[1]]]),
                         category=names(unique_peaks),
                         fill=unname(unlist(conditions_colour_code)),
                         col=NA,
                         #cat.pos=c(0,0),
                         cat.dist = c(0,0))
  grid.arrange(gTree(children=g), top="Binding Site Overlaps", bottom=gsub('/', '', opt$result_dir))
  invisible(capture.output(dev.off()))
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
  
  png(paste(output_prefix, 'raw_consensus_peaks.png', sep=""))
  plt = plot(euler(e), main=gsub('/', '', opt$result_dir), quantities=TRUE, fills=unname(unlist(conditions_colour_code)))
  invisible(capture.output(ggsave(filename=paste(output_prefix, 'raw_consensus_peaks.png', sep=''), plot=plt)))
  invisible(capture.output(dev.off()))
}


png(paste(output_prefix, 'raw_consensus_heatmap.png', sep=""))
dba.plotHeatmap(dbObj.consensus)
invisible(capture.output(dev.off()))

png(paste(output_prefix, 'raw_pca_condition.png', sep=""))
dba.plotPCA(dbObj, masks=!dbObj.total$masks$Consensus, attributes=DBA_CONDITION, label=DBA_ID, vColors=(colours))
invisible(capture.output(dev.off()))
if(length(unique(dbObj$samples$Factor)) > 1){
  png(paste(output_prefix, 'raw_pca_factor.png', sep=""))
  dba.plotPCA(dbObj, masks=!dbObj.total$masks$Consensus, attributes=DBA_FACTOR, label=DBA_ID)
  invisible(capture.output(dev.off()))
}

# Plot peaks over genome
plt <- covplot(c(unique_peaks, shared_peaks), title="Peaks over Genome") + 
  scale_color_manual(values=rev(c(unlist(unname(conditions_colour_code[1:length(unique_peaks)])), 'grey'))) + 
  scale_fill_manual(values=rev(c(unlist(unname(conditions_colour_code[1:length(unique_peaks)])), 'grey')))
invisible(capture.output(ggsave(filename=paste(output_prefix, 'genome_peaks.png', sep=''), plot=plt, dpi=320)))
plt <- plt + facet_grid(chr ~ .id)
invisible(capture.output(ggsave(filename=paste(output_prefix, 'genome_peaks_split.png', sep=''), plot=plt, dpi=320)))

# Plot peaks related to TSS sites
tagMatrices <- list()
for (p in names(unique_peaks)){
  tagMatrix <- getTagMatrix(unique_peaks[[p]], windows=promoters)
  if (length(tagMatrix) == 0){
    cat("No peaks at promoter sites for", p, "\n")
    break
  }else{
    tagMatrices[[p]] <- tagMatrix
    cat(dim(tagMatrix)[[1]], "peaks at promoter sites for", p, "\n")
    plt <- tagHeatmap(tagMatrix, 
               xlab="bp at TSS", 
               ylab="Peaks", 
               title=paste(dim(tagMatrix)[[1]],'Peaks at Promoters', p, sep=" - "),
               palette=if_else(conditions_colour_code[[p]] == "#00BFC4", 'Greens', 'Reds'), 
    )
    invisible(capture.output(ggsave(paste(output_prefix, 'raw_TSS_heatmap_', p, '_peaks.png', sep=''), plot=plt, dpi=320)))
  }
}

# Plot TSS profile of peaks
plt <- plotAvgProf(tagMatrices, xlim=c(-3000, 3000), conf=0.95, resample=1000) +
  scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))) +
  scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)])))
invisible(capture.output(ggsave(paste(output_prefix, 'raw_TSS_profile_peaks.png', sep=''), plot=plt, dpi=320)))

peaks <- c(unique_peaks, shared_peaks)

# Output to bed files
for (p in names(peaks)){
  write.table(as.data.frame(peaks[[p]]), file=paste(output_prefix, p, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
}

# ========= Get Annotations =========
peakAnnoList <- list()
for (p in names(peaks)){
  cat("\nAnnotating", p)
  anno <- annotatePeak(peaks[[p]], 
                       TxDb=txdb,
                       annoDb=annoDb)
  
  plt <- upsetplot(anno, vennpie=TRUE) + ggtitle(p)
  invisible(capture.output(ggsave(paste(output_prefix, '_', p, '_annotated_peaks.png', sep=''), plot=plt, dpi=320, bg='white')))
  
  peakAnnoList[[p]] <- anno
  
  # Write annotation to file
  write.table(anno, file=paste(output_prefix, 'annotated_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
  
  tryCatch(
    {
      genes <- list()
      genes[[p]] <- anno@anno$geneId
      names(genes) = sub("_", "\n", names(genes))
      
      # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
      compKEGG <- compareCluster(geneCluster=genes,
                                 fun="enrichKEGG",
                                 pvalueCutoff=0.05,
                                 pAdjustMethod="BH",
                                 organism=keggOrg) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
      if (!is.null(compKEGG)){
        #plt <- dotplot(compKEGG, showCategory = 10, title = "KEGG Pathway Enrichment Analysis")
        plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=conditions_colour_code[[p]], n=15)
        invisible(capture.output(ggsave(filename=paste(output_prefix, p, '_annotated_kegg_analysis.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
        
        # Write annotations to csv
        write.table(as.data.frame(compKEGG), file=paste(output_prefix, p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
      }
    },error = function(e)
    {
      message(e)
    }
  )
  
  for (ont in c('ALL', 'CC', 'MF', 'BP')){
    tryCatch(
      {
        genes <- list()
        genes[[p]] <- anno@anno$SYMBOL
        names(genes) = sub("_", "\n", names(genes))
        
        compGO <- compareCluster(geneCluster=genes,
                                 keyType='SYMBOL',
                                 OrgDb=annoDb,
                                 fun="enrichGO",
                                 ont=ont,
                                 pvalueCutoff=0.05,
                                 pAdjustMethod="BH",
                                 readable=TRUE
        ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
        if (!is.null(compGO)){
          compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
          #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
          plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=conditions_colour_code[[p]], n=15)
          invisible(capture.output(ggsave(filename=paste(output_prefix, '_', p, '_annotated_go_analysis', '_', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
          
          # Write annotations to csv
          write.table(as.data.frame(compGO), file=paste(output_prefix, '_', p, 'annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
        }
      },error = function(e)
      {
        message(e)
      }
    )
    
  }
}

plt <- plotAnnoBar(peakAnnoList)
invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_distribution.png', sep=''), plot=plt, dpi=320)))

plt <- plotDistToTSS(peakAnnoList)
invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_TSS_distribution.png', sep=''), plot=plt, dpi=320)))


# ========= END OF OCCUPANCY ANALYSIS =========




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
output_prefix <- gsub('.csv', '_', paste(paste(result_dir, sep=''), samplesheet, sep=""))
output_prefix <- gsub("diffbind_samplesheet_", "", output_prefix)

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
  
  dbObj.counted <- dba.count(dbObj.caller_consensus, bUseSummarizeOverlaps=TRUE, 
                             peaks=consensus_peaks, 
                             minOverlap=1, 
                             score=DBA_SCORE_NORMALIZED,
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
  dbObj.counted <- dba.count(dbObj.noblacklist, bUseSummarizeOverlaps=TRUE, 
                             peaks=consensus_peaks, 
                             minOverlap=1, 
                             score=DBA_SCORE_NORMALIZED,
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
    # Normalize ("safest method")
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
png(paste(output_prefix, 'consensus_peaks_counted_normalized_heatmap.png', sep=''))
dba.plotHeatmap(dbObj.norm)
invisible(capture.output(dev.off()))
png(paste(output_prefix, 'consensus_peaks_counted_normalized_pca.png', sep=''))
dba.plotPCA(dbObj.norm, attributes=DBA_CONDITION, label=DBA_ID, vColors=(colours))
invisible(capture.output(dev.off()))

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
dbObj.analyzed <- dba.analyze(dbObj.contrast, method=DBA_ALL_METHODS)
cat("After analyzing:\n")
dbObj.analyzed


# Plots
tryCatch(
  {
    png(paste(output_prefix, 'analyzed_venn.png', sep=''))
    dba.plotVenn(main="DE Binding Sites Identified by Method",
                 dbObj.analyzed, 
                 contrast=1, 
                 method=DBA_ALL_METHODS
    )
    invisible(capture.output(dev.off()))
  }, error=function(e){
    message("No figure\n", e)
    invisible(capture.output(dev.off()))
    invisible(file.remove(paste(output_prefix, 'analyzed_venn.png', sep='')))
  }
)

for (m in c(DBA_DESEQ2, DBA_EDGER)){
  
  tryCatch(
    {
      #png(paste(output_prefix, 'analyzed_venn_', m, '.png', sep=''))
      #dba.plotVenn(dbObj.analyzed, method=m, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
      #invisible(capture.output(dev.off()))
      tryCatch(
        {
          v <- dba.plotVenn(dbObj.analyzed, method=m, contrast=1, 
                            bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE,
                            main=paste("Binding Sites - ", m, sep=''))
        }, error=function(e){
          message('Peaksets do not meet specified criteria for venn')
        }
      )
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
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      invisible(file.remove(paste(output_prefix, 'analyzed_venn_', m, '.png', sep='')))
    }
  )
  tryCatch(
    {
      png(paste(output_prefix, 'analyzed_volcano_', m, '.png', sep=''))
      par(mfrow=c(1,2))
      dba.plotVolcano(dbObj.analyzed, method=m)
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      invisible(file.remove(paste(output_prefix, 'analyzed_volcano_', m, '.png', sep='')))
    }
  )
  tryCatch(
    {
      png(paste(output_prefix, 'analyzed_ma_', m, '.png', sep=''))
      dba.plotMA(dbObj.analyzed, method=m)
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      invisible(file.remove(paste(output_prefix, 'analyzed_ma_', m, '.png', sep='')))
    }
  )
  tryCatch(
    {
      png(paste(output_prefix, 'analyzed_box_', m, '.png', sep=''))
      dba.plotBox(dbObj.analyzed, method=m)
      invisible(capture.output(dev.off()))
    }, error=function(e){
      message("No figure\n", e)
      invisible(capture.output(dev.off()))
      invisible(file.remove(paste(output_prefix, 'analyzed_box_', m, '.png', sep='')))
    }
  )
}

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
    invisible(file.remove(paste(output_prefix, 'analyzed_profile.png', sep='')))
  }
)

# ========= Generate DE Analysis Reports =========
reports <- list()
reports[["DESeq2"]] <- dba.report(dbObj.analyzed, method=DBA_DESEQ2, contrast=1, th=1)
reports[["edgeR"]] <- dba.report(dbObj.analyzed, method=DBA_EDGER, contrast=1, th=1)

for (report in names(reports)){
  # Report columns are seqnames, start, end, width, strand, Conc, Conc_Group1, Conc_Group2, Fold, p.value, FDR
  # Create bed files for each keeping only significant peaks (p<0.05)
  # Comparing those whereby Fold > 0 vs Fold < 0, indicating enrichment gain vs loss of group 1 over group 2
  
  out <- as.data.frame(reports[[report]])
  #out <- out[as.data.frame(findOverlaps(both, reports[[report]]))[["subjectHits"]], ]
  
  cat("\n", report, "Report:\n")
  print(head(out))
  print(tail(out))
  
  
  gained <- out %>% 
    dplyr::filter(FDR < 0.05 & Fold > 0) %>% 
    dplyr::select(seqnames, start, end)
  lost <- out %>% 
    dplyr::filter(FDR < 0.05 & Fold < 0) %>% 
    dplyr::select(seqnames, start, end)
  
  cat('\n', dim(out)[[1]], 'peaks in', report, 'report:\n')
  cat('\n\t', dim(gained)[[1]] + dim(lost)[[1]], 'statistically significant peaks', '( FDR <', dbObj.analyzed$config$th, ')\n')
  cat('\n\t', dim(gained)[[1]], 'DE peaks in', dbObj.contrast$contrasts[[1]]$name1, '(Fold-change > 0)', '\n')
  cat('\n\t', dim(lost)[[1]], 'DE peaks in', dbObj.contrast$contrasts[[1]]$name2, '(Fold-change < 0)', '\n')
  
  # Write complete report to file
  write.table(out, file=paste(output_prefix, 'analyzed_report_', report, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
  
  # Write DE result to bed files
  write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  
  
  # Plot profile heatmaps for all significant (FDR < 0.05) sites for each method
  tryCatch(
    {
      profile_colors <- list()
      for (i in 1:length(unique(dbObj.analyzed$samples$Condition))) {
        profile_colors[[unique(dbObj.analyzed$samples$Condition)[i]]] <- c('white', colours[[i]])
      }
      #profile_colors <- rev(profile_colors)
      # Plots all significant sites among both conditions (will include signals)
      profiles_significant <- dba.plotProfile(dbObj.analyzed, merge=c(DBA_REPLICATE), normalize=TRUE, 
                                              sites=GRanges(out %>% dplyr::filter(FDR < 0.05))
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
      invisible(file.remove(paste(output_prefix, 'analyzed_profiles_', report,'.png', sep='')))
    }
  )
  
  gained <- out %>% 
    dplyr::filter(FDR < 0.05 & Fold > 0)
  lost <- out %>% 
    dplyr::filter(FDR < 0.05 & Fold < 0)
  
  # Write complete DE result to files
  write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  
  # Use absolute fold-change for plotting magnitude
  gained$Fold <- abs(gained$Fold)
  lost$Fold <- abs(lost$Fold)
  
  gpeaks <- GenomicRanges::GRangesList(Gained=gained, Lost=lost)
  names(gpeaks) <- c(dbObj.contrast$contrasts[[1]]$name1, dbObj.contrast$contrasts[[1]]$name2)
  
  # Plot peaks gained/lost over genome between conditions
  tryCatch(
      {
        plt <- covplot(gpeaks, title=paste("Peaks over Genome", report, sep=' - '), weightCol='Fold') + 
          scale_color_manual(values=rev(c(colours[1:length(unique_peaks)]))) + 
          scale_fill_manual(values=rev(c(colours[1:length(unique_peaks)])))
        invisible(capture.output(ggsave(filename=paste(output_prefix, '_', report, '_significant_merged_peaks.png', sep=''), plot=plt, dpi=320)))
        plt <- plt + facet_grid(chr ~ .id)
        invisible(capture.output(ggsave(filename=paste(output_prefix, '_', report, '_significant_peaks.png', sep=''), plot=plt, dpi=320)))
        }, error=function(e){
          message("No figure\n", e)
      }
    )
  
  # ========= Get Annotations =========
  peakAnnoList <- list()
  for (p in names(gpeaks)){
    if (length(gpeaks[[p]]) > 0){
      cat("\nAnnotating", p, "\n")
      anno <- annotatePeak(gpeaks[[p]], 
                           TxDb=txdb,
                           annoDb=annoDb)
      peakAnnoList[[p]] <- anno
      
      plt <- upsetplot(anno, vennpie=TRUE) + ggtitle(p)
      invisible(capture.output(ggsave(paste(output_prefix, '_', report, '_', p, '_annotated_peaks.png', sep=''), plot=plt, dpi=320, bg='white')))
      
      # Write annotation to file
      write.table(anno, file=paste(output_prefix, 'annotated_', report, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
  
      tryCatch(
        {
          genes <- list()
          genes[[p]] <- anno@anno$geneId
          names(genes) = sub("_", "\n", names(genes))
          
          # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
          compKEGG <- compareCluster(geneCluster=genes,
                                     fun="enrichKEGG",
                                     pvalueCutoff=0.05,
                                     pAdjustMethod="BH",
                                     organism=keggOrg) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          
          if (!is.null(compKEGG)){
            plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=conditions_colour_code[[p]], n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_kegg_analysis.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
            
            # Write annotations to csv
            write.table(as.data.frame(compKEGG), file=paste(output_prefix, report, '_', p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
          }
        },error = function(e)
        {
          message(e)
        }
      )
      
      for (ont in c('ALL', 'CC', 'MF', 'BP')){
        tryCatch(
          {
            compGO <- compareCluster(geneCluster=genes,
                                     OrgDb=annoDb,
                                     fun="enrichGO",
                                     ont=ont,
                                     pvalueCutoff=0.05,
                                     pAdjustMethod="BH",
                                     readable=TRUE
            ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
            
            if (!is.null(compGO)){
              compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
              #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
              plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=conditions_colour_code[[p]], n=15)
              invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_go_analysis', '_', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
              
              # Write annotations to csv
              write.table(as.data.frame(compGO), file=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
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
  
  
  plt <- plotAnnoBar(peakAnnoList)
  invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))
  
  plt <- plotDistToTSS(peakAnnoList)
  invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_TSS_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))
  
}


# For considering sites identified by DESeq2 AND edgeR
#reports[["DESeq2_and_edgeR"]] <- c(reports[["DESeq2"]], reports[["edgeR"]])
#names(reports[["DESeq2_and_edgeR"]]) <- 1:length(names(reports[["DESeq2_and_edgeR"]]))
df_both <- merge(as.data.frame(reports[["DESeq2"]]), as.data.frame(reports[["edgeR"]]), by=c("seqnames", "start", "end"), suffixes=c(".DESeq2", ".edgeR"))
reports[["DESeq2_and_edgeR"]] <- GRanges(df_both)
report <- "DESeq2_and_edgeR"

#both <- dba.plotVenn(dbObj.analyzed, contrast=1, method=DBA_ALL_METHODS, bDB=TRUE)$inAll
#out_both <- out[as.data.frame(findOverlaps(both, reports[["DESeq2_and_edgeR"]]))[["subjectHits"]], ]

out <- as.data.frame(reports[[report]])
#out <- out[as.data.frame(findOverlaps(both, reports[[report]]))[["subjectHits"]], ]

cat("\n", report, "Report:\n")
print(head(out))
print(tail(out))


gained <- out %>% 
  dplyr::filter(FDR.DESeq2 < 0.05 & Fold.DESeq2 > 0 & FDR.edgeR < 0.05 & Fold.edgeR > 0) %>% 
  dplyr::select(seqnames, start, end)
lost <- out %>% 
  dplyr::filter(FDR.DESeq2 < 0.05 & Fold.DESeq2 < 0 & FDR.edgeR < 0.05 & Fold.edgeR < 0) %>% 
  dplyr::select(seqnames, start, end)

cat('\n', dim(out)[[1]], 'peaks in', report, 'report:\n')
cat('\n\t', dim(gained)[[1]] + dim(lost)[[1]], 'statistically significant peaks', '( FDR <', dbObj.analyzed$config$th, ')\n')
cat('\n\t', dim(gained)[[1]], 'DE peaks (positive fold-change) in', dbObj.contrast$contrasts[[1]]$name1, '\n')
cat('\n\t', dim(lost)[[1]], 'DE peaks (negative fold-change) in', dbObj.contrast$contrasts[[1]]$name2, '\n')

# Write complete report to file
write.table(out, file=paste(output_prefix, 'analyzed_report_', report, '.tsv', sep=''), sep="\t", quote=F, row.names=F)

# Write DE result to bed files
write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)

gained <- out %>% 
  dplyr::filter(FDR.DESeq2 < 0.05 & Fold.DESeq2 > 0 & FDR.edgeR < 0.05 & Fold.edgeR > 0)
lost <- out %>% 
  dplyr::filter(FDR.DESeq2 < 0.05 & Fold.DESeq2 < 0 & FDR.edgeR < 0.05 & Fold.edgeR < 0)

# Write complete DE result to files
write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
# Use absolute fold-change for plotting magnitude
#gained$Fold <- abs(gained$Fold)
#lost$Fold <- abs(lost$Fold)

gpeaks <- GenomicRanges::GRangesList(Gained=gained, Lost=lost)
names(gpeaks) <- c(dbObj.contrast$contrasts[[1]]$name1, dbObj.contrast$contrasts[[1]]$name2)

# Plot peaks gained/lost over genome between conditions
tryCatch(
  {
    plt <- covplot(gpeaks, title=paste("Peaks over Genome", report, sep=' - '), ) + #weightCol='Fold') + 
      scale_color_manual(values=rev(c(colours[1:length(unique_peaks)]))) + 
      scale_fill_manual(values=rev(c(colours[1:length(unique_peaks)])))
    invisible(capture.output(ggsave(filename=paste(output_prefix, '_', report, '_significant_merged_peaks.png', sep=''), plot=plt, dpi=320)))
    plt <- plt + facet_grid(chr ~ .id)
    invisible(capture.output(ggsave(filename=paste(output_prefix, '_', report, '_significant_peaks.png', sep=''), plot=plt, dpi=320)))
  }, error=function(e){
    message("No figure\n", e)
    invisible(capture.output(dev.off()))
  }
)


# ========= Get Annotations =========
peakAnnoList <- list()
for (p in names(gpeaks)){
  if (length(gpeaks[[p]]) > 0){
    cat("\nAnnotating", p, "\n")
    anno <- annotatePeak(gpeaks[[p]], 
                         TxDb=txdb,
                         annoDb=annoDb)
    peakAnnoList[[p]] <- anno
    
    plt <- upsetplot(anno, vennpie=TRUE) + ggtitle(p)
    invisible(capture.output(ggsave(paste(output_prefix, '_', report, '_', p, '_annotated_peaks.png', sep=''), plot=plt, dpi=320, bg='white')))
    
    # Write annotation to file
    write.table(anno, file=paste(output_prefix, 'annotated_', report, '_', p, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
    
    tryCatch(
      {
        genes <- list()
        genes[[p]] <- anno@anno$geneId
        names(genes) = sub("_", "\n", names(genes))
        
        # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
        compKEGG <- compareCluster(geneCluster=genes,
                                   fun="enrichKEGG",
                                   pvalueCutoff=0.05,
                                   pAdjustMethod="BH",
                                   organism=keggOrg) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
        
        if (!is.null(compKEGG)){
          plt <- make_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=conditions_colour_code[[p]], n=15)
          invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_kegg_analysis.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
          
          # Write annotations to csv
          write.table(as.data.frame(compKEGG), file=paste(output_prefix, report, '_', p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
        }
      },error = function(e)
      {
        message(e)
      }
    )
    
    for (ont in c('ALL', 'CC', 'MF', 'BP')){
      tryCatch(
        {
          compGO <- compareCluster(geneCluster=genes,
                                   OrgDb=annoDb,
                                   fun="enrichGO",
                                   ont=ont,
                                   pvalueCutoff=0.05,
                                   pAdjustMethod="BH",
                                   readable=TRUE
          ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          
          if (!is.null(compGO)){
            compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
            #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
            plt <- make_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=conditions_colour_code[[p]], n=15)
            invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_go_analysis', '_', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
            
            # Write annotations to csv
            write.table(as.data.frame(compGO), file=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
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

plt <- plotAnnoBar(peakAnnoList)
invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))

plt <- plotDistToTSS(peakAnnoList)
invisible(capture.output(ggsave(filename=paste(output_prefix, 'peaks_annotation_TSS_distribution_', report, '.png', sep=''), plot=plt, dpi=320)))

cat('\nFINISHED!\n')
