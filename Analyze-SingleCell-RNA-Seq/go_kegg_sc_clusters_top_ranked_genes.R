library(clusterProfiler)
library(paletteer)
library(stringr)
library(ggplot2)
library(optparse)

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

make_anno_dotplot <- function(df, title="", ylabel="Description", colour="#4393C3", n=15,
                              title_size=16,
                              text_size=2,
                              axis_title_size=8,
                              axis_x_size=10,
                              axis_y_size=10,
                              #legend_key_size=1.6,
                              legend_title_size=14,
                              legend_text_size=12,
                              colour_text=FALSE){
  # Expects df input from clusterProfiler output as dataframe
  
  # Example usage:
  # plt <- make_anno_dotplot(compGO@compareClusterResult, title="YourTitle", ylabel="GO Term", n=15)
  # ggsave(filename='YourFigure.png', plot=plt, dpi=400, units='mm', width=200, height=200)
  
  df <- head(df[order(df$p.adjust, df$FoldEnrichment, -xtfrm(df$GeneRatio), -xtfrm(df$BgRatio), df$Description), ], n=n)
  df$Order <- 1:dim(df)[1]
  df$ycolour <- "black"
  if ("ONTOLOGY" %in% colnames(df)){
    if (length(unique(df$ONTOLOGY)) > 1){
      df$Description <- paste(df$ONTOLOGY, df$Description, sep=' - ')
      if (colour_text){
        df$ycolour <- ifelse(grepl("BP -", df$Description), "#4393C3", df$ycolour)
        df$ycolour <- ifelse(grepl("CC -", df$Description), "#D6604D", df$ycolour)
        df$ycolour <- ifelse(grepl("MF -", df$Description), "#639863FF", df$ycolour)
      }
    }
  }
  
  # Re-format y-axis labels to not squish graph
  for (d in 1:length(df$Description)){
    i <- 1
    s <- df$Description[d]
    s <- str_remove(s, " - Mus musculus \\(house mouse\\)")
    df$Description[d] <- ""
    while (i < length(strsplit(s, ' ')[[1]]) + 1){
      df$Description[d] <- paste(df$Description[d], paste(strsplit(s, ' ')[[1]][i:(i+3)], collapse = ' '), sep='\n')
      i <- (i+4)
    }
    df$Description[d] <- gsub(" NA", "", df$Description[d])
    df$Description[d] <- substring(df$Description[d], 2, nchar(df$Description[d]))
  }
  #numGR <- sapply(df$GeneRatio, function(x) eval(parse(text = x)))
  # Plot
  plt <- ggplot() +
    geom_point(data=head(df, n=n),
               aes(x = -log(p.adjust), 
                   y = reorder(Description, -Order),# -p.adjust), 
                   colour = Count, 
                   size = unname(unlist(sapply(GeneRatio, function(x) eval(parse(text=x)))))*100,
               ),
    ) + 
    theme_classic() +
    theme(#axis.text.y = element_text(colour=rev(df$ycolour), face = "bold", size = axis_y_size),
      axis.text.y = element_text(colour = "black", face = "bold", size = axis_y_size),
      axis.title.x = element_text(size = axis_x_size*1.4),
      axis.text.x = element_text(size = axis_x_size, colour='black'),
      axis.title.y = element_text(size = axis_y_size*1.4),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
      legend.title = element_text(size = legend_title_size, face = "bold"),
      legend.text = element_text(size = legend_text_size)
    ) +
    scale_color_gradient(low = "black", high = colour) +
    ggtitle(title)
  plt$labels$x <- "-log(p.adjust)"
  plt$labels$y <- ylabel
  plt$labels$size <- "GenePercentage"
  plt$labels$colour <- "GeneCount"
  return(plt)
}

fig_size <- function(type = c("single", "onehalf", "double", "full"),
                     aspect = 0.75,  # height/width
                     pointsize = 8) {
  type <- match.arg(type)
  w <- switch(type,
              single  = 3.5,
              onehalf = 5.0,
              double  = 7.0,
              full = 9.5)
  h <- w * aspect
  list(width = w, height = h, pointsize = pointsize)
}

option_list = list(
  make_option(c("-i", "--input"), type="character", default="ranked_genes.csv", help=".csv file from scanpy ranked genes groups columns as GENE_NAME_0, PVAL_0, GENE_NAME_1, PVAL_1, etc...for each cluster", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default="hg38", help="Assembly to annotate genes/peaks (e.g. hg19, hg38, mm9, mm10, rn6)", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="PEA_Results/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-d", "--database"), type="character", default="ucsc", help="Database reference for peaks gene annotations, ucsc (default) or ensembl", metavar="character"),
  make_option(c("-l", "--annotation_level"), type="character", default="transcript", help="Level parameter for annotatePeak, 'gene' or 'transcript'", metavar="character"),
  make_option(c("--lfc"), type="double", default=0.585, help="Magnitude of log2foldchange to define significant up/down regulation of genes", metavar="double"),
  make_option(c("--padj"), type="double", default=0.05, help="Significance threshold (p.adjust value) for DEGs", metavar="double"),
  make_option(c("--david_user"), type="character", default="earezza@ohri.ca", help="User email for DAVID web tools (must be registered, https://david.ncifcrf.gov/content.jsp?file=DAVID_WebService.html)", metavar="character"),
  make_option(c("--minGSSize"), type="integer", default=10, help="minimal size of genes annotated for testing", metavar="integer"),
  make_option(c("--maxGSSize"), type="integer", default=500, help="maximal size of genes annotated for testing", metavar="integer"),
  make_option(c("--colours_discrete"), type="character", default="khroma::muted", help="Palette from paletteer for discrete colours, see https://pmassicotte.github.io/paletteer_gallery/#discrete-palettes", metavar="character"),
  make_option(c("--reverse_d_palette"), type="logical", action="store_true", default=FALSE, help="Reverse the discrete palette colours", metavar="logical"),
  make_option(c("--colours_continuous"), type="character", default="ggthemes::Classic Blue", help="Palette from paletteer for continuous colours, see https://pmassicotte.github.io/paletteer_gallery/#continuous-palettes", metavar="character"),
  make_option(c("--reverse_c_palette"), type="logical", action="store_true", default=FALSE, help="Reverse the continuous palette colours", metavar="logical"),
  make_option(c("-c", "--colours"), type="character", default="clusters_hexcolours.csv", help=".csv file defining cluster hex colours with columns as CLUSTER, HEX", metavar="character"),
  make_option(c("--figsize"), type="character", default="full", help="Sizing for figures, options are 'full' for full-page width, 'single' for single-column width, 'onehalf' for 1.5 column width, 'double' for double-column width", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (str_sub(opt$result_dir, -1) != "/"){
  opt$result_dir = cat(opt$result_dir, '/', sep='')
}
if (!file.exists(opt$input)){
  cat(opt$input, "does not exist...check path and current working directory.")
  q()
}
if (!(opt$assembly %in% c('mm10', 'mm9', 'hg38', 'hg19', 'rn6'))){
  cat(opt$assembly, "not a valid choice. Only supports mm9, mm10, hg19, hg38, rn6 assemblies.")
  q()
}

figs <- fig_size(opt$figsize, aspect = 0.85, pointsize = 10)

# Set discrete colour palette
colours_discrete <- tryCatch(
  {
    # Try user‑specified palette
    paletteer::paletteer_d(opt$colours_discrete, direction=ifelse(opt$reverse_d_palette, -1, 1))
  },
  error = function(e) {
    message(
      "Error occurred: ", conditionMessage(e),
      "\nDefaulting to 'khroma::muted' discrete colour palette\n"
    )
    # Fallback palette
    paletteer::paletteer_d("khroma::muted", direction=ifelse(opt$reverse_d_palette, -1, 1))
  }
)
# Add colours to palette if not enough
if (length(colours_discrete) < 20){
  muted_plus_9 <- c(
    "#7C7C8C",  # muted blue‑grey
    "#B48A76",  # muted warm tan
    "#5F6A3A",  # dark moss green
    "#6C8A9C",  # muted steel blue
    "#A27D9C",  # dusty mauve
    "#8C5E4A",  # muted umber
    "#6A7F5C",  # sage‑olive
    "#A3A0B3",  # soft lilac‑grey
    "#C3A39A"   # pale clay
  )
  class(muted_plus_9) <- class(colours_discrete)
  colours_discrete <- c(colours_discrete, muted_plus_9[1:(20-length(colours_discrete))])
  class(colours_discrete) <- class(muted_plus_9)
}

colours_continuous <- tryCatch(
  {
    # Try user‑specified palette
    paletteer::paletteer_c(opt$colours_continuous, n=100, direction=ifelse(opt$reverse_c_palette, -1, 1))
  },
  error = function(e) {
    message(
      "Error occurred: ", conditionMessage(e),
      "\nDefaulting to 'ggthemes::Classic Blue' continuous colour palette\n"
    )
    # Fallback palette
    paletteer::paletteer_c("ggthemes::Classic Blue", n=100, direction=ifelse(opt$reverse_c_palette, -1, 1))
  }
)


setwd(dirname(opt$input))
# Set directories
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}
cat("Output files will be in", opt$result_dir, "\n")

# Get assembly annotation info
anno_ref <- load_annotation(opt$assembly, opt$database)

# Read in differential expression ranked genes list
df <- read.csv(opt$input, header = TRUE, check.names = FALSE)
#colnames(df)

# Extract number of clusters and their respective gene lists
gene_cols <- colnames(df)[grepl('_n', colnames(df))]
p_cols <- colnames(df)[grepl('_pvals_adj', colnames(df))]
l_cols <- colnames(df)[grepl('_logfoldchanges', colnames(df))]

if (length(p_cols) == 0){
  cat("No pvals_adj columns...using pvals...")
  p_cols <- colnames(df)[grepl('_pvals', colnames(df))]
}
# sanity check
if (length(gene_cols) != length(p_cols)){
  cat("Number of clusters mismatch for columns of gene names and p values...quitting.")
  q()
}

cat("\nObtaining significant genes for each cluster...\n")
cluster_genes <- list()
gsea_cluster_genes <- list()
for (c in 1:length(gene_cols)){
  cat('\nCluster', c - 1)
  df_cluster <- df[c(gene_cols[c], p_cols[c], l_cols[c])]
  
  cat("\nGenes before significance filter:", dim(df_cluster)[1])
  df_cluster <- df_cluster[df_cluster[colnames(df_cluster)[2]] <= opt$padj, ]
  cat("\nGenes after significance filter:", dim(df_cluster)[1])
  cluster_genes[[as.character(c-1)]] <- df_cluster[[colnames(df_cluster)[1]]]
  
  geneList <- df[l_cols[c]][[l_cols[c]]]
  names(geneList) <- df[gene_cols[c]][[gene_cols[c]]]
  geneList <- sort(geneList, decreasing = TRUE)
  gsea_cluster_genes[[as.character(c-1)]] <- geneList

}

cluster_colours <- list()
if (is.null(opt$colours)){
  for (c in names(cluster_genes)){
    cluster_colours[[c]] <- colours_discrete[as.integer(c)]
  }
} else {
  hex <- read.csv(opt$colours, header = TRUE, check.names = FALSE)
  for (c in names(cluster_genes)){
    cluster_colours[[c]] <- hex[ hex[, 1] == c, ][[colnames(hex)[2]]]
  }
}

if (!is.null(opt$gene_universe)){
  bg_genes <- read.csv(opt$gene_universe, header=FALSE)$V1
  
  if (opt$assembly == "hg19" | opt$assembly == "hg38"){
    bg_entrez <- mapIds(org.Hs.eg.db, keys = bg_genes, column = "ENTREZID", keytype = "SYMBOL")
  }else if (opt$assembly == "mm9" | opt$assembly == "mm10"){
    bg_entrez <- mapIds(org.Mm.eg.db, keys = bg_genes, column = "ENTREZID", keytype = "SYMBOL")
  }else if (opt$assembly == "rn6"){
    bg_entrez <- mapIds(org.Rn.eg.db, keys = bg_genes, column = "ENTREZID", keytype = "SYMBOL")
  }
}


# Annotate significant gene sets for each cluster
tryCatch(
  {
    # # Output log to file
    # con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    # sink(con, split = FALSE)                 # normal output
    # #sink(con, type = "message", split = TRUE)  # messages
    cat("\n\n---------- Performing GO + KEGG ----------\n")
    
    pdf(paste(opt$result_dir, 'GO_KEGG_clusters_top_differential_genes.pdf', sep=""),
        width  = figs$width,
        height = figs$height,
        pointsize = figs$pointsize)  #8 or 9–10)
    
    for (c in names(cluster_genes)){
      cat("\n === Cluster", c, "===\n")
      # Set plot colors
      colour <- cluster_colours[[c]]
      heat_colour <- cluster_colours[[c]]
      
      # Mapper for EntrezID to gene SYMBOL
      df_mapper <- data.frame()
      genes_entrez <- list()
      
      if (opt$assembly == "hg19" | opt$assembly == "hg38"){
        genes_entrez[[c]] <- mapIds(org.Hs.eg.db, keys = cluster_genes[[c]], column = "ENTREZID", keytype = "SYMBOL")
      }else if (opt$assembly == "mm9" | opt$assembly == "mm10"){
        genes_entrez[[n]] <- mapIds(org.Mm.eg.db, keys = cluster_genes[[c]], column = "ENTREZID", keytype = "SYMBOL")
      }else if (opt$assembly == "rn6"){
        genes_entrez[[n]] <- mapIds(org.Rn.eg.db, keys = cluster_genes[[c]], column = "ENTREZID", keytype = "SYMBOL")
      }
      df_mapper <- rbind(df_mapper, data.frame(geneId=unname(genes_entrez[[c]]), SYMBOL=cluster_genes[[c]]))
      
      # Mapper for EntrezID to gene SYMBOL
      mapper <- df_mapper[!duplicated(df_mapper), ]
      
      tryCatch(
        {
          genes <- list()
          genes[[c]] <- mapper$geneId
          names(genes) = sub("_", "\n", names(genes))
          
          # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
          compKEGG <- compareCluster(geneCluster=genes,
                                     fun="enrichKEGG",
                                     pvalueCutoff=opt$padj,
                                     pAdjustMethod="BH",
                                     minGSSize = opt$minGSSize,
                                     maxGSSize = opt$maxGSSize,
                                     organism=anno_ref$keggOrg,
                                     universe = bg_entrez
          ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          if (!class(compKEGG) == 'compareClusterResult'){
            cat("\nNo KEGG results.\n")
            #next
          }else{
            cat('\n', dim(compKEGG@compareClusterResult)[1], 'KEGG results\n')
          }
          if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
            # Map EntrezIDs to gene SYMBOL
            compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
            myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
            for (i in 1:length(myEntrez)){
              compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
            }
            # Write annotations to csv
            df <- as.data.frame(compKEGG@compareClusterResult)
            df <- df[order(df$p.adjust, df$FoldEnrichment, -xtfrm(df$GeneRatio), -xtfrm(df$BgRatio), df$Description), ]
            write.table(df, file=paste(opt$result_dir, c, '_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
            
            plt <- make_anno_dotplot(compKEGG@compareClusterResult, 
                                     title=paste('KEGG - Cluster ', c, sep=""), 
                                     ylabel="KEGG Category", 
                                     colour=cluster_colours[[c]], 
                                     n=15,
                                     title_size=figs$pointsize*1.6,
                                     text_size=figs$pointsize/3,
                                     axis_title_size=figs$pointsize*0.8,
                                     axis_x_size=figs$pointsize*1.2,
                                     axis_y_size=figs$pointsize*1.2,
                                     #legend_key_size=figs$pointsize/5,
                                     legend_title_size=figs$pointsize*1.4,
                                     legend_text_size=figs$pointsize*1.2
            )
            print(plt)
            #invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], p, '_consensus_annotated_KEGG.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
            
          } else{
            cat("\nNo annotation results\n")
            remove(compKEGG)
            gc()
          }
        },error = function(e)
        {
          message(e)
          remove(compKEGG)
          gc()
        }
      )
      invisible(capture.output(gc()))
      
      
      #for (ont in c('ALL', 'CC', 'MF', 'BP')){
      for (ont in c('ALL', 'CC', 'MF', 'BP')){
        tryCatch(
          {
            genes <- list()
            genes[[c]] <- mapper$SYMBOL
            names(genes) = sub("_", "\n", names(genes))
            #cat("\nGetting ", p, ' GO ', ont, '\n')
            
            compGO <- compareCluster(geneCluster=genes,
                                     keyType='SYMBOL',
                                     OrgDb=anno_ref$annoDb,
                                     fun="enrichGO",
                                     ont=ont,
                                     pvalueCutoff=opt$padj,
                                     pAdjustMethod="BH",
                                     minGSSize = opt$minGSSize,
                                     maxGSSize = opt$maxGSSize,
                                     readable=TRUE,
                                     universe = bg_genes
            ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
            if (!class(compGO) == 'compareClusterResult'){
              cat("\nNo GO results", "for", ont, ".\n")
              next
            }else{
              cat('\n', dim(compGO@compareClusterResult)[1], 'GO', ont, 'results\n')
            }
            if ((!is.null(compGO)) & (dim(compGO@compareClusterResult)[1] > 0)){
              #compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology # found instances of inaccuracies...
              # Write annotations to csv
              df <- as.data.frame(compGO@compareClusterResult)
              df <- df[order(df$p.adjust, df$FoldEnrichment, -xtfrm(df$GeneRatio), -xtfrm(df$BgRatio), df$Description), ]
              write.table(df, file=paste(opt$result_dir, c, '_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_anno_dotplot(compGO@compareClusterResult, 
                                       title=paste("GO (", ont, ") - Cluster ", c, sep=""), 
                                       ylabel="GO Term", 
                                       colour=cluster_colours[[c]], 
                                       n=15,
                                       title_size=figs$pointsize*1.6,
                                       text_size=figs$pointsize/3,
                                       axis_title_size=figs$pointsize*0.8,
                                       axis_x_size=figs$pointsize*1.2,
                                       axis_y_size=figs$pointsize*1.2,
                                       #legend_key_size=figs$pointsize/5,
                                       legend_title_size=figs$pointsize*1.4,
                                       legend_text_size=figs$pointsize*1.2
              )
              #invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], p, '_annotated_GO-', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
              print(plt)
              remove(compGO)
              gc()
            } else{
              cat("\nNo GO", ont, "annotation results\n")
              remove(compGO)
              gc()
            }
          },error = function(e)
          {
            message(e)
            remove(compGO)
            gc()
          }
        )
        invisible(capture.output(gc()))
      }
      
      # GSEA
      for (ont in c('ALL', 'CC', 'MF', 'BP')){
        cat("\nGSEA - ", ont, "\n")
        tryCatch(
          {
            gsea <- gseGO(geneList=na.omit(gsea_cluster_genes[[c]]), 
                          ont = ont, 
                          keyType = "SYMBOL", 
                          minGSSize = opt$minGSSize,
                          maxGSSize = opt$maxGSSize, 
                          pvalueCutoff = opt$pvalue, 
                          verbose = TRUE, 
                          OrgDb = anno_ref$annoDb, 
                          pAdjustMethod = "BH",
                          nPermSimple = 10000,
                          eps = 0
            )
            if (!class(gsea) == 'gseaResult'){
              cat("\nNo GSEA -", ont, "results.\n")
              next
            }else{
              cat('\n', dim(gsea@result)[1], 'GSEA', ont, 'results\n')
            }
            if ((!is.null(gsea)) & (dim(gsea@result)[1] > 0)){
              gsea@result <- gsea@result[order(gsea@result$p.adjust, decreasing=FALSE),] # Sort by most signiicant
              gsea@result$SYMBOL <- gsea@result$core_enrichment
              # Write annotations to csv
              df <- as.data.frame(gsea$result)
              write.table(df, file=paste(opt$result_dir, c, '_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_anno_dotplot(gsea$result, 
                                       title=paste("GSEA (", ont, ") - Cluster ", c, sep=""), 
                                       ylabel="Term", 
                                       colour=cluster_colours[[c]], 
                                       n=15,
                                       title_size=figs$pointsize*1.6,
                                       text_size=figs$pointsize/3,
                                       axis_title_size=figs$pointsize*0.8,
                                       axis_x_size=figs$pointsize*1.2,
                                       axis_y_size=figs$pointsize*1.2,
                                       #legend_key_size=figs$pointsize/5,
                                       legend_title_size=figs$pointsize*1.4,
                                       legend_text_size=figs$pointsize*1.2
              )
              remove(gsea)
              gc()
            } else{
              cat("\nNo annotation results\n")
              remove(gsea)
              gc()
            }
            
          },error = function(e)
          {
            message(e)
            remove(gsea)
            gc()
          }
        ) 
      }
      
      
      # genes_entrez <- list()
      # genes_entrez[[p]] <- anno@anno$geneId
      # names(genes_entrez) = sub("_", "\n", names(genes_entrez))
      # for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
      #   cat("\nGetting ", p, ' DAVID ', annotation_type, ' annotations...\n')  
      #   # DAVID Annotation
      #   tryCatch(
      #     {
      #       cat("\nDAVID - ", annotation_type, "\n")
      #       entrez <- unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))])
      #       
      #       py_require(c("pandas"))
      #       py_require(c("suds"))
      #       py_run_string("import pandas as pd")
      #       py_run_string("import sys")
      #       py_run_string("from suds.client import Client")
      #       
      #       # create a service client using the wsdl.
      #       py_run_string("client = Client('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl')")
      #       py_run_string("client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')")
      #       
      #       #authenticate user email
      #       py_run_string("client.service.authenticate(r.opt['david_user'])")
      #       
      #       # Read input gene list file, convert ids to a comma-delimited string and upload the list to DAVID
      #       py_run_string("inputIds = ','.join(r.entrez)")
      #       py_run_string("client.service.addList(inputIds, 'ENTREZ_GENE_ID', r.p, 0)")
      #       
      #       # setCategories
      #       py_run_string("categorySting = str(client.service.setCategories(r.annotation_type))")
      #       
      #       #getChartReport
      #       py_run_string("thd = r.opt['fdr']")
      #       py_run_string("ct = r.opt['minGSSize']")
      #       py_run_string("chartReport = client.service.getChartReport(thd,ct)")
      #       py_run_string("chartRow = len(chartReport)")
      #       py_run_string("print ('Total chart records:',chartRow)")
      #       
      #       if (py$chartRow > 0){
      #         if (annotation_type == "KEGG_PATHWAY"){
      #           splitter <- ":"
      #         }else{
      #           splitter <- "~"
      #         }
      #         # parse chartReport
      #         records <- data.frame(
      #           ID            = character(),
      #           Category      = character(),
      #           Description   = character(),
      #           GeneRatio     = character(),
      #           BgRatio       = character(),
      #           pvalue        = numeric(),
      #           p.adjust      = numeric(),
      #           FDR           = numeric(),
      #           geneID        = character(),
      #           Count         = integer(),
      #           foldEnrichment= numeric(),
      #           id            = character(),
      #           stringsAsFactors = FALSE
      #         )
      #         
      #         for (simpleChartRecord in py$chartReport) {
      #           df <- data.frame(
      #             ID            = strsplit(simpleChartRecord$termName, splitter)[[1]][1],
      #             Category      = simpleChartRecord$categoryName,
      #             Description   = strsplit(simpleChartRecord$termName, splitter)[[1]][2],
      #             GeneRatio     = paste0(simpleChartRecord$listHits, "/", simpleChartRecord$listTotals),
      #             BgRatio       = paste0(simpleChartRecord$popHits, "/", simpleChartRecord$popTotals),
      #             pvalue        = simpleChartRecord$ease,
      #             p.adjust      = simpleChartRecord$benjamini,
      #             FDR           = simpleChartRecord$afdr,
      #             geneID        = gsub(", ", "/", simpleChartRecord$geneIds),
      #             Count         = simpleChartRecord$listHits,
      #             foldEnrichment= simpleChartRecord$foldEnrichment,
      #             id            = simpleChartRecord$id,
      #             stringsAsFactors = FALSE
      #           )
      #           records <- rbind(records, df)
      #         }
      #         # py_run_string("records = pd.DataFrame()")
      #         # py_run_string("for simpleChartRecord in chartReport:
      #         #                 df = pd.DataFrame(index=[records.shape[0]], data={
      #         #                     'ID' : simpleChartRecord['termName'].split(r.splitter)[0],
      #         #                     'Category' : simpleChartRecord['categoryName'],
      #         #                     'Description' : simpleChartRecord['termName'].split(r.splitter)[1],
      #         #                     'GeneRatio': str(simpleChartRecord['listHits']) + '/' + str(simpleChartRecord['listTotals']), 
      #         #                     'BgRatio': str(simpleChartRecord['popHits']) + '/' + str(simpleChartRecord['popTotals']), 
      #         #                     'pvalue' : simpleChartRecord['ease'],
      #         #                     'p.adjust' : simpleChartRecord['benjamini'],
      #         #                     'FDR' : simpleChartRecord['afdr'],
      #         #                     'geneID' : simpleChartRecord['geneIds'].replace(', ', '/'),
      #         #                     'Count' : simpleChartRecord['listHits'],
      #         #                     'foldEnrichment' : simpleChartRecord['foldEnrichment'],
      #         #                     'id' : simpleChartRecord['id']
      #         #                     })
      #         #                 records = pd.concat([records, df])")
      #         # #py_run_string("records.reset_index(drop=True, inplace=True)")
      #         #compD <- py$records
      #         compD <- records
      #         colnames(compD)[colnames(compD) == "foldEnrichment"] <- "FoldEnrichment"
      #         qobj <- tryCatch(qvalue(p=as.numeric(compD$pvalue), lambda=opt$pvalue, pi0.method="bootstrap"),
      #                          error=function(e) NULL)
      #         if (class(qobj) == "qvalue") {
      #           qvalues <- qobj$qvalues
      #         } else {
      #           qvalues <- NA
      #         }
      #         compD$qvalue <- qvalues
      #         
      #         compDAVID <- new("enrichResult",
      #                          result         = compD,
      #                          pvalueCutoff   = 1,
      #                          pAdjustMethod  = "BH",
      #                          organism       = opt$assembly,
      #                          ontology       = annotation_type,
      #                          gene           = entrez,
      #                          keytype        = "ENTREZ_GENE_ID")
      #         rm(compD)
      #         
      #         if (!class(compDAVID) == 'enrichResult'){
      #           cat("\nNo DAVID", annotation_type, "results.\n")
      #           next
      #         }else{
      #           cat('\n', dim(compDAVID@result)[1], 'DAVID', annotation_type,  'results\n')
      #         }
      #         # old deprecated function
      #         #compDAVID <- enrichDAVID(
      #         #  unname(genes_entrez[[n]][!is.na(unname(genes_entrez[[n]]))]),
      #         #  idType = "ENTREZ_GENE_ID",
      #         #  minGSSize = opt$minGSSize,
      #         #  maxGSSize = opt$maxGSSize,
      #         #  annotation = annotation_type,
      #         #  pvalueCutoff = opt$pvalue,
      #         #  pAdjustMethod = "BH",
      #         #  #species = NA,
      #         #  david.user=opt$david_user
      #         #)
      #         
      #         if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
      #           # Map EntrezIDs to gene SYMBOL
      #           compDAVID@result$SYMBOL <- compDAVID@result$geneID
      #           myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
      #           for (i in 1:length(myEntrez)){
      #             compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
      #           }
      #           # Write annotations to csv
      #           df <- as.data.frame(compDAVID@result)
      #           df <- df[order(df$p.adjust, df$FoldEnrichment, -xtfrm(df$GeneRatio), -xtfrm(df$BgRatio), df$Description), ]
      #           write.table(df, file=paste(result_dirs[[p]], p, '_consensus_annotated_DAVID_', annotation_type, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
      #           
      #           plt <- make_anno_dotplot(compDAVID@result, 
      #                                    title=paste('DAVID - ', p, " (Consensus)", sep=""), 
      #                                    ylabel=paste(annotation_type,"Category", sep=' '), 
      #                                    colour=colours_continuous, 
      #                                    n=15,
      #                                    title_size=figs$pointsize*1.6,
      #                                    text_size=figs$pointsize/3,
      #                                    axis_title_size=figs$pointsize*0.8,
      #                                    axis_x_size=figs$pointsize*1.2,
      #                                    axis_y_size=figs$pointsize*1.2,
      #                                    #legend_key_size=figs$pointsize/5,
      #                                    legend_title_size=figs$pointsize*1.4,
      #                                    legend_text_size=figs$pointsize*1.2
      #           )
      #           #invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], 'DAVID_annotation_', annotation_type, '_', p, '_dotplot.png', sep=''), plot=plt, dpi=320)))
      #           print(plt)
      #           remove(compDAVID)
      #           gc()
      #         } else{
      #           cat("\nNo DAVID", annotation_type, "annotation results\n")
      #           remove(compDAVID)
      #           gc()
      #         }
      #       }
      #     },error = function(e)
      #     {
      #       message(e)
      #       remove(compDAVID)
      #       gc()
      #     }
      #   )
      # }
      
    }
 
    # # Upset plot for annotated genes
    # upsetlist <- list()
    # for (p in names(peakAnnoList)) {
    #   upsetlist[[p]] <- peakAnnoList[[p]]@anno$SYMBOL
    # }
    # upset_colors <- list()
    # for (n in names(peakAnnoList)){
    #   upset_colors[[ conditions_colour_code[[n]] ]] <- length(unique(peakAnnoList[[n]]@anno$SYMBOL))
    # }
    # upset_colors <- sort(unlist(upset_colors), decreasing=TRUE)
    # 
    # plt <- upset(fromList(upsetlist), 
    #              order.by = "freq", 
    #              nsets = length(names(peakAnnoList)),
    #              sets.bar.color = names(upset_colors),
    #              empty.intersections = "on",
    #              set_size.show = TRUE,
    #              set_size.angles = 0,
    #              set_size.scale_max = dim(fromList(upsetlist))[[1]]*1.25,
    #              sets.x.label = "Gene Set Size",
    #              mainbar.y.label = "Intersection Size of Gene Sets",
    #              mb.ratio = c(0.65, 0.35),
    #              #main.bar.color = unlist(conditions_colour_code),
    #              text.scale = c(figs$pointsize/6, figs$pointsize/6, figs$pointsize/7, figs$pointsize/7, figs$pointsize/6, figs$pointsize/6)
    #              # order: intersection size title, set size title,
    #              # intersection tick labels, set tick labels,
    #              # intersection bar labels, set size labels
    # )
    # print(plt)
    # grid.text("Annotated Genes in Consensus Peaksets",x = 0.65, y=0.95, gp=gpar(fontsize=figs$pointsize*1.2))
    
    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    # #sink()                  # stop messages
    # sink()                  # stop normal output
    # close(con)
  }
  
)

