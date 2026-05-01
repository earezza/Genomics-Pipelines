# -----------------------------------------------------------
# Script Name: Analyze Peaks
# Purpose: This script performs differential peaks analysis using DiffBind.
# Author: Eric Arezza
# Date: 2024-04-23
# Version: 1.0
# -----------------------------------------------------------

# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
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
  library(paletteer)
  library(pheatmap)
  library(grid)
  library(reticulate)
  library(qvalue)
  library(optparse)
}))

# Define Functions
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
  
  df <- head(df[order(df$p.adjust, -xtfrm(df$GeneRatio), -xtfrm(df$BgRatio), df$Description), ], n=n)
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

make_anno_piebar <- function(df, 
                             type='pie', 
                             title="Distribution of Sites", 
                             specific=TRUE, 
                             colours=c(paletteer_d("khroma::muted" ), "#7C7C8CFF", "#B48A76FF"),
                             title_size=10,
                             text_size=1.8,
                             axis_title_size=8,
                             axis_x_size=5,
                             axis_y_size=5,
                             legend_key_size=1.6,
                             legend_title_size=6,
                             legend_text_size=6
                             ){
  # Expects df input from annotatePeak anno@anno as dataframe ( or list of dataframes for stacked bar )
  
  # Example usage:
  # plt <- make_anno_piebar(as.data.frame(anno@anno), type='pie', specific=TRUE, title="YourTitle")
  # ggsave(filename="YourFigure.png", plot=plt, dpi=400, units='mm', width=170, height=130)
  if (specific){
    df$Group <- ifelse(grepl("Promoter \\(<=1kb\\)", df$annotation), 'Promoter\n(<=1kb)', 'Non-promoter')
    df$Group <- ifelse(grepl("Promoter \\(1-2kb\\)", df$annotation), 'Promoter\n(1-2kb)', df$Group)
    df$Group <- ifelse(grepl("Promoter \\(2-3kb\\)", df$annotation), 'Promoter\n(2-3kb)', df$Group)
  } else {
    df$Group <- ifelse(grepl("Promoter", df$annotation), 'Promoter', 'Non-promoter')
  }
  
  df$Group <- ifelse(grepl("Exon", df$annotation), 'Exon', df$Group)
  df$Group <- ifelse(grepl("Intron", df$annotation), 'Intron', df$Group)
  df$Group <- ifelse(grepl("5' UTR", df$annotation), "5' UTR", df$Group)
  df$Group <- ifelse(grepl("3' UTR", df$annotation), "3' UTR", df$Group)
  df$Group <- ifelse(grepl("Downstream", df$annotation), "Downstream\n(<=300bp)", df$Group)
  df$Group <- ifelse(grepl("Distal Intergenic", df$annotation), "Distal\nIntergenic", df$Group)
  # Note: "Enhancer" typically not output from annotatePeak...
  df$Group <- ifelse(grepl("Enhancer", df$annotation), 'Enhancer', df$Group)
  
  colour_map <- c(
    "3' UTR" = colours[1],
    "5' UTR" = colours[2],
    "Distal\nIntergenic" = colours[3],
    "Downstream\n(<=300bp)" = colours[4],
    "Exon" = colours[5],
    "Intron" = colours[6],
    "Promoter\n(<=1kb)" = colours[7],
    "Promoter\n(1-2kb)" = colours[8],
    "Promoter\n(2-3kb)" = colours[9],
    "Promoter" = colours[7]
  )
  colour_map <- colour_map[names(colour_map) %in% unique(df$Group)]
  #df$Group <- factor(df$Group, levels = names(colour_map))
  
  # Data to plot
  df_plot <- df %>% count(Group)
  df_plot$Frequency <- round( 100 * (df_plot$n / sum(df_plot$n)), 2)
  df_plot <- df_plot[order(df_plot$Frequency, decreasing=TRUE), ]
  rownames(df_plot) <- 1:nrow(df_plot)
  colnames(df_plot) <- c('Region', 'Count', 'Frequency')
  df_plot$Region <- as.factor(df_plot$Region)
  
  # Colours
  df_plot$Colour <- unname(colour_map[df_plot$Region])
  #myColors <- colours[1:length(unique(df_plot$Region))]
  #names(myColors) <- levels(df_plot$Region)
  custom_colors_fill <- scale_fill_manual(values = colour_map,#myColors, 
                                          name = "Region",
                                          labels = str_c(df_plot$Region, ' (', df_plot$Frequency, '%)', sep=''))
  
  # Plot Bar Chart
  if (type == 'bar'){
    plt <- ggplot(df_plot, aes(x=reorder(Region, -Count), y=Count, fill=Region)) + 
      geom_bar(stat='identity') +
      custom_colors_fill +
      geom_text(stat='identity', aes(label=Count), vjust=-1, size = text_size) +
      theme_classic() +
      theme(axis.text.x = element_text(size = axis_x_size, color = "black")) +
      theme(axis.text.y = element_text(size = axis_y_size, color = "black")) +
      theme(axis.title = element_text(size = axis_title_size, color = "black")) +
      scale_y_continuous(limits = c(0, max(df_plot$Count) * 1.1)) +
      guides(fill="none") +
      ggtitle(title) + 
      theme(plot.title = element_text(size = title_size))
    plt$labels$x <- "Region"
    plt$labels$y <- "Number of Sites"
  }
  else if (type == 'pie'){
    plt <- ggplot(df_plot, aes(x='', y=Frequency, fill=reorder(Region, -Frequency))) +
      geom_bar(stat = "identity", position = "stack") +
      coord_polar(theta = "y", start=0) +
      custom_colors_fill +
      theme_void() +
      guides(fill='legend') +
      ggtitle(title) + 
      theme(plot.title = element_text(size = title_size)) +
      theme(plot.background = element_rect(color = 'white', fill = "white")) + 
      theme(plot.margin = margin(0,1,0,0, "cm")) +
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(
        legend.key.size = unit(legend_key_size, "in"),
        legend.title = element_text(size = legend_title_size),
        legend.text = element_text(size = legend_text_size)
      )
  }
  return(plt)
}

# make_anno_stackedbar <- function(df_list, title="Distribution of Sites", specific=TRUE, colours=paletteer_d("khroma::muted"), xlabel="Percentage", ylabel="Group"){
#   # Expects list of named dataframes where each is an output from annotatePeak (anno@anno as dataframe)
#   
#   # plt <- make_anno_stackedbar(anno_list, specific=TRUE, title="YourTitle", ylabel="SampleGroup")
#   # ggsave(filename="YourFigure.png", plot=plt, dpi=400, units='mm', width=180, height=120)
#   
#   
#   for (n in names(df_list)){
#     df <- df_list[[n]]
#     if (specific){
#       df$Group <- ifelse(grepl("Promoter \\(<=1kb\\)", df$annotation), 'Promoter\n(<=1kb)', 'Non-promoter')
#       df$Group <- ifelse(grepl("Promoter \\(1-2kb\\)", df$annotation), 'Promoter\n(1-2kb)', df$Group)
#       df$Group <- ifelse(grepl("Promoter \\(2-3kb\\)", df$annotation), 'Promoter\n(2-3kb)', df$Group)
#     } else {
#       df$Group <- ifelse(grepl("Promoter", df$annotation), 'Promoter', 'Non-promoter')
#     }
#     
#     df$Group <- ifelse(grepl("Exon", df$annotation), 'Exon', df$Group)
#     df$Group <- ifelse(grepl("Intron", df$annotation), 'Intron', df$Group)
#     df$Group <- ifelse(grepl("5' UTR", df$annotation), "5' UTR", df$Group)
#     df$Group <- ifelse(grepl("3' UTR", df$annotation), "3' UTR", df$Group)
#     df$Group <- ifelse(grepl("Downstream", df$annotation), "Downstream\n(<=300bp)", df$Group)
#     df$Group <- ifelse(grepl("Distal Intergenic", df$annotation), "Distal\nIntergenic", df$Group)
#     # Note: "Enhancer" typically not output from annotatePeak...
#     df$Group <- ifelse(grepl("Enhancer", df$annotation), 'Enhancer', df$Group)
#     
#     # Data to plot
#     df_plot <- df %>% count(Group)
#     df_plot$Frequency <- round( 100 * (df_plot$n / sum(df_plot$n)), 2)
#     df_plot <- df_plot[order(df_plot$Frequency, decreasing=TRUE), ]
#     rownames(df_plot) <- 1:nrow(df_plot)
#     colnames(df_plot) <- c('Region', 'Count', 'Frequency')
#     df_plot$Region <- as.factor(df_plot$Region)
#     df_list[[n]] <- df_plot
#   }
#   
#   my_plots <- df_list
#   combined_df <- bind_rows(my_plots, .id = "group")
#   # Order from highest to lowest
#   total_freqs <- combined_df %>% 
#     group_by(Region) %>% 
#     summarise(across(c("Frequency"), sum, na.rm = TRUE))
#   total_freqs <- total_freqs[order(total_freqs$Frequency, decreasing=FALSE), ]
#   
#   # Colours
#   myColors <- colours[1:length(total_freqs$Region)]
#   #names(myColors) <- fct_infreq(levels(combined_df$Region))
#   custom_colors_fill <- scale_fill_manual(values = myColors, 
#                                           name = "Region",
#                                           labels = total_freqs$Region
#   )
#   # Plot
#   plt <- ggplot(combined_df, aes(x = fct_rev(group), y = Frequency, fill = factor(Region, levels = total_freqs$Region) )) +
#     geom_bar(stat = "identity", position="stack") +
#     custom_colors_fill +
#     theme_classic() +
#     theme(axis.text.x = element_text(size = 12, color = "black")) +
#     theme(axis.text.y = element_text(size = 12, color = "black")) +
#     theme(axis.title = element_text(size = 14, color = "black")) +
#     ggtitle(title) + 
#     theme(plot.title = element_text(size = 16, color = "black")) + 
#     coord_flip() +
#     labs(x = ylabel, y = xlabel, fill = "Region") +
#     theme(plot.title = element_text(hjust = 0.5)) + 
#     theme(
#       legend.key.size = unit(0.75, 'cm'),
#       legend.title = element_text(size = 12),
#       legend.text = element_text(size = 10)
#     ) + 
#     guides(fill = guide_legend(reverse = TRUE))
#   
#   return(plt)
# }

# make_dotplot <- function(df, 
#                          title="", 
#                          ylabel="Description", 
#                          colour="#56B1F7", 
#                          n=15,
#                          title_size,
#                          axis_title_size,
#                          axis_x_size,
#                          axis_y_size,
#                          legend_title_size,
#                          legend_text_size){
#   df$ycolour <- "black"
#   if ("ONTOLOGY" %in% colnames(df)){
#     df$Description <- paste(df$ONTOLOGY, df$Description, sep=' - ')
#     df$ycolour <- ifelse(grepl("BP -", df$Description), 'blue', df$ycolour)
#     df$ycolour <- ifelse(grepl("CC -", df$Description), 'red', df$ycolour)
#     df$ycolour <- ifelse(grepl("MF -", df$Description), 'darkgreen', df$ycolour)
#   }
#   df <- df[order(df$p.adjust, decreasing=FALSE),]
#   
#   # Re-format y-axis labels to not squish graph
#   for (d in 1:length(df$Description)){
#     i <- 1
#     s <- df$Description[d]
#     s <- str_remove(s, " - Mus musculus \\(house mouse\\)")
#     df$Description[d] <- ""
#     while (i < length(strsplit(s, ' ')[[1]]) + 1){
#       df$Description[d] <- paste(df$Description[d], paste(strsplit(s, ' ')[[1]][i:(i+5)], collapse = ' '), sep='\n')
#       i <- (i+6)
#     }
#     df$Description[d] <- gsub(" NA", "", df$Description[d])
#     df$Description[d] <- substring(df$Description[d], 2, nchar(df$Description[d]))
#   }
#   
#   # Plot
#   plt <- ggplot() +
#     geom_point(data=head(df, n=n),
#                aes(x = -log(p.adjust), 
#                    y = reorder(Description, -p.adjust), 
#                    colour = Count, 
#                    size = unname(unlist(sapply(GeneRatio, function(x) eval(parse(text=x)))))*100,
#                ),
#     ) + 
#     theme_classic() +
#     theme(axis.text.y = element_text(colour=rev(head(df$ycolour, n=n)))) +
#     scale_color_gradient(low = "black", high = colour) +
#     ggtitle(title) 
#   plt$labels$x <- "-log(p.adjust)"
#   plt$labels$y <- ylabel
#   plt$labels$size <- "GenePercentage"
#   plt$labels$colour <- "GeneCount"
#   return(plt)
# }

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

# ======= Get command-line optional arguments =======
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, help="DiffBind-formatted sample sheet", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default="mm10", help="Assembly to annotate genes/peaks (e.g. hg19, hg38, mm9, mm10, rn6)", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="Peaks_Analysis/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-d", "--database"), type="character", default="ucsc", help="Database reference for peaks gene annotations, ucsc (default) or ensembl", metavar="character"),
  make_option(c("-l", "--annotation_level"), type="character", default="transcript", help="Level parameter for annotatePeak, 'gene' or 'transcript'", metavar="character"),
  make_option(c("--combine_callers"), type="logical", action="store_true", default=FALSE, help="Flag to add peaks from callers instead of taking consensus peaks", metavar="logical"),
  make_option(c("--combine_replicates"), type="logical", action="store_true", default=FALSE, help="Flag to add peaks from all replicates instead of taking consensus peaks", metavar="logical"),
  make_option(c("-b", "--blacklisted_keep"), type="logical", action="store_true", default=FALSE, help="Flag to keep blacklisted regions in raw peaks files", metavar="logical"),
  make_option(c("--lfc"), type="double", default=0.585, help="Magnitude of log2foldchange to define significant up/down enrichment of binding sites", metavar="double"),
  make_option(c("--fdr"), type="double", default=0.05, help="Significance threshold (false discovery rate, a.k.a. p.adjust value)", metavar="double"),
  make_option(c("--occupancy_only"), type="logical", action="store_true", default=FALSE, help="Flag to only perform peaks occupancy analysis", metavar="logical"),
  make_option(c("--david_user"), type="character", default="earezza@ohri.ca", help="User email for DAVID web tools (must be registered, https://david.ncifcrf.gov/content.jsp?file=DAVID_WebService.html)", metavar="character"),
  make_option(c("--minGSSize"), type="integer", default=10, help="minimal size of genes annotated for testing", metavar="integer"),
  make_option(c("--maxGSSize"), type="integer", default=500, help="maximal size of genes annotated for testing", metavar="integer"),
  make_option(c("--bg_genes"), type="character", default=NULL, help="Custom .txt file of 1 column list all genes to use as background set in pathway analysis", metavar="character"),
  make_option(c("--figsize"), type="character", default="full", help="Sizing for figures, options are 'full' for full-page width, 'single' for single-column width, 'onehalf' for 1.5 column width, 'double' for double-column width", metavar="character"),
  make_option(c("--colours_discrete"), type="character", default="khroma::muted", help="Palette from paletteer for discrete colours, see https://pmassicotte.github.io/paletteer_gallery/#discrete-palettes", metavar="character"),
  make_option(c("--reverse_d_palette"), type="logical", action="store_true", default=FALSE, help="Reverse the discrete palette colours", metavar="logical"),
  make_option(c("--colours_continuous"), type="character", default="ggthemes::Classic Blue", help="Palette from paletteer for continuous colours, see https://pmassicotte.github.io/paletteer_gallery/#continuous-palettes", metavar="character"),
  make_option(c("--reverse_c_palette"), type="logical", action="store_true", default=FALSE, help="Reverse the continuous palette colours", metavar="logical"),
  make_option(c("-p", "--python"), type="character", default="base", help="Name of python environment (can be conda or virtualenv) to use for running DAVID annotations", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Load available Python environment
tryCatch(
  {
    use_condaenv(opt$python, required = TRUE)
    print("Conda environment found.\n")
  },
  error = function(e) {
    warning(e)
  }
)
tryCatch(
  {
    use_virtualenv(opt$python, required = TRUE)
    print("Virtualenv found.\n")
  },
    error = function(e) {
      warning(e)
    }
)
py_config()

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
if (length(colours_discrete) < 11){
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
  colours_discrete <- c(colours_discrete, muted_plus_9[1:(11-length(colours_discrete))])
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

if (!is.null(opt$bg_genes)){
  bg_genes <- read.csv(opt$bg_genes, header=FALSE)$V1
  
  if (opt$assembly == "hg19" | opt$assembly == "hg38"){
    bg_entrez <- mapIds(org.Hs.eg.db, keys = bg_genes, column = "ENTREZID", keytype = "SYMBOL")
  }else if (opt$assembly == "mm9" | opt$assembly == "mm10"){
    bg_entrez <- mapIds(org.Mm.eg.db, keys = bg_genes, column = "ENTREZID", keytype = "SYMBOL")
  }else if (opt$assembly == "rn6"){
    bg_entrez <- mapIds(org.Rn.eg.db, keys = bg_genes, column = "ENTREZID", keytype = "SYMBOL")
  }
} else {
  bg_genes <- NULL
  bg_entrez <- NULL
}

# Setup directories
if (file.exists(opt$file)) {
  setwd(dirname(opt$file))
} else {
  cat('\nCannot find DiffBind samplesheet --file', opt$file, "\n")
  quit()
}
if (str_sub(opt$result_dir, -1) != "/"){
  opt$result_dir = cat(opt$result_dir, '/', sep='')
}
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}
result_dir <- paste(opt$result_dir, "Occupancy_Analysis/", sep='')
if (!file.exists(result_dir)) {
  dir.create(result_dir)
}


tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "wt")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    
    cat("Run options:\n")
    for (i in which(names(opt) != "help")) {
      cat(names(opt)[i], '=', paste(opt)[i], "\n")
    }
    cat("\nlog2FC of", opt$lfc, "equates to FC of", round(2^0.585, 2), '\n')
    
    if (!(opt$assembly %in% c('mm10', 'mm9', 'hg38', 'hg19', 'rn6'))){
      cat(opt$assembly, "not a valid choice. Only supports mm9, mm10, hg19, hg38, rn6 assemblies.")
      sink()                  # stop normal output
      close(con)
      quit()
    }
    
    cat(
    "\n\n=============================== START OF OCCUPANCY ANALYSIS ===============================
    Here, peaks declared by peak caller(s) are used to identify differential binding (DB) 
    between conditions. Overlapping peaks between peak sets (.bed files) are determined by 
    the range of the peaks where 1bp overlaps are merged into a consensus peak range. 
    Blacklisted regions (due to known genomic regions and sequencing artifacts) are by
    default removed from the peaks, but can be kept using the --blacklisted_keep flag.
    
    By default, a final consensus peak set is formed whereby consensus peaks found in at 
    least 2 replicates and 2/3 peak callers form the condition's final consensus peak set.
    Including all replicate and caller peaks can be done with the --combine_replicates and 
    --combine_callers flags.
    
    Occupancy analysis alone provides a general analysis for DB by simply considering and 
    comparing where peaks exist. 
    
    Following this (if chosen to execute), affinity analysis can provide more conservative 
    and statistically validated results for DB since the enrichment of read counts are 
    accounted for (peak shapes and signal strength) among the identified consensus regions
    from the occupancy analysis.\n\n"
    )
    
    
    # ========= SETUP RUN AND VARIABLES =========
    # ========= Load peaksets =========
    samplesheet <- basename(opt$file)
    sample_table <- read.csv(samplesheet)
    #output_prefix <- gsub('.csv', '_', paste(result_dir, samplesheet, sep=""))
    #output_prefix <- gsub("diffbind_samplesheet_", "", output_prefix)
    
    # Get annotations reference and respective promoter regions
    anno_ref <- load_annotation(opt$assembly, opt$database)
    promoters <- getPromoters(TxDb=anno_ref$txdb, upstream=3000, downstream=3000)
    
    # Variables to get average fragment size for each sample, used later for counts
    fragment_size <- 1:length(read.csv(samplesheet)$SampleID)
    frag_sizes <- list()
    sample_reps <- unique(paste0(sample_table$Condition, ':', sample_table$Replicate))
    sample_reps_split <- str_split(sample_reps, pattern = ':')
    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
)


# Load peaks, filter blacklisted regions, setup variables, plot some QC
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    
    pdf(paste(result_dir, 'SupplementaryQC.pdf', sep=""),
        width  = figs$width,
        height = figs$height,
        pointsize = figs$pointsize)  # or 9–10)
    
    for (s in 1:length(sample_reps)){
      mean_fragment_size <- average_fragment_length(unique(sample_table[sample_table$Condition == sample_reps_split[[s]][1] 
                                                                        & sample_table$Replicate == sample_reps_split[[s]][2], ]$bamReads)
      )
      title(main = sample_reps[s])
      frag_sizes[[sample_reps[s]]] = mean_fragment_size
      fragment_size[as.integer(rownames(sample_table[sample_table$Condition == sample_reps_split[[s]][1] 
                                                     & sample_table$Replicate == sample_reps_split[[s]][2], ]))] <- mean_fragment_size
    }
    #invisible(capture.output(dev.off())) 
    #invisible(capture.output(gc())) 
    
    #for (b in unique(read.csv(samplesheet)$Condition)){
    #  for (r in unique(read.csv(samplesheet)$Replicate)){
    #    png(paste(result_dir, 'fragment_length_', b, '-', r, '.png', sep=""))
    #    mean_fragment_size <- average_fragment_length(read.csv(samplesheet)$bamReads[[which(read.csv(samplesheet)$Condition == b)[1]]], plot=TRUE)
    #    for (i in which(read.csv(samplesheet)$Condition == b & read.csv(samplesheet)$Replicate == r)){
    #      fragment_size[i] <- mean_fragment_size
    #      frag_sizes[[b]] <- mean_fragment_size
    #    }
    #    invisible(capture.output(dev.off())) 
    #  }
    #}
    #invisible(capture.output(gc())) 
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
    cat("\n\n---------- Raw peaksets ----------\n\n")
    print(dbObj)
    
    # Colour codes for consistency in plots (add more colours if needed)
    # colours <- c("#00BFC4", "#F8766D", "#7CAE00", "#C77CFF", "#e69e02", "#00A9FF", "#C77CFF", "#FF61CC", 
    #              "#FF0000", "#FF4D00", "#80FF00", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2" , "#D55E00" , "#CC79A7"
    # )
    
    conditions_colour_code <- list()
    for (i in 1:length(unique(dbObj$samples$Condition))) {
      conditions_colour_code[[unique(dbObj$samples$Condition)[i]]] <- colours_discrete[i]
    }
    conditions_colour_code[['Shared']] <- "grey"
    
    #png(filename=paste(result_dir, 'raw_heatmap.png', sep=''))
    dba.plotHeatmap(dbObj, margin=15, cexRow = 0.8, cexCol = 0.8)
    #title(main=c("Correlation Heatmap - Raw Peaks"))
    mtext(c("Correlation Heatmap - Raw Peaks"), side = 1, line = 2)
    #invisible(capture.output(dev.off()))
    #invisible(capture.output(gc()))
    
    # Show overlap rates for each condition
    cat("\nTotal peaks overlapped when found in at least (1, 2, ...) replicates/callers for each condition:\n")
    #png(paste(result_dir, "raw_overlap_rates.png", sep=""))
    par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
    #par(mfrow=c(length(unique(dbObj$samples$Condition)), 1), mar = c(5, 5, 4, 4))
    for (c in unique(dbObj$samples$Condition)) {
      cat('\n', c, '\n')
      olap.rate <- dba.overlap(dbObj, dba.mask(dbObj, attribute=DBA_CONDITION, value=c, combine='or'), mode=DBA_OLAP_RATE)
      #cat(olap.rate)
      #cat('\n')
      plot(olap.rate, type='l', xlim=c(0.75, length(olap.rate)), ylab='# overlapping peaks', xlab='# peaksets (replicates and peak callers)', col=conditions_colour_code[[c]])
      text(x=1:length(olap.rate), y=olap.rate, olap.rate, cex=0.75)
      axis(side=1, at=1:length(olap.rate))
      title(main=c("Overlap Rate - Raw Peaks", c))
    }
    #invisible(capture.output(dev.off()))
    #invisible(capture.output(gc()))
    
    # ========= Remove Blacklisted Regions =========
    # Remove blacklisted regions to ignore irrelevant peaks (blacklisted regions from ENCODE, genome selected is based on prediction from bam files)
    tryCatch (
      {
        if (!opt$blacklisted_keep){
          cat("\n\n---------- Removing blacklisted regions ----------\n\n")
          sink(con, append = TRUE, type = "message")
          dbObj.noblacklist <- dba.blacklist(dbObj, blacklist=TRUE, greylist=FALSE)
          sink(type = "message")
          
          #blacklisted_peaks <- dba.blacklist(dbObj.noblacklist, Retrieve=DBA_BLACKLISTED_PEAKS)
          cat("\n\n---------- After blacklisted regions removed ----------\n\n")
          print(dbObj.noblacklist)
        }else{
          cat("\n\n---------- Blacklisted regions not removed, proceeding with raw peaksets...----------\n\n")
          dbObj.noblacklist <- dbObj
        }
      },error = function(e)
      {
        message(e)
      }
    )
    # if (!exists("dbObj.noblacklist")) {
    #   cat("\nBlacklisted regions not removed, proceeding with raw peaksets...\n")
    #   dbObj.noblacklist <- dbObj
    #   print(dbObj.noblacklist)
    # }
    
    #png(paste(result_dir, 'raw_noblacklist_heatmap.png', sep=''))
    dba.plotHeatmap(dbObj.noblacklist, margin=15, cexRow = 0.8, cexCol = 0.8)
    #title(main=c("Correlation Heatmap - Blacklist Regions Removed"))
    mtext(c("Correlation Heatmap - Blacklist Regions Removed"), side = 1, line = 2)
    #invisible(capture.output( dev.off() ))
    #invisible(capture.output(gc()))
    
    # Show overlap rates for each condition
    cat("\nTotal peaks overlapped when found in at least (1, 2, ...) replicates/callers for each condition:\n")
    #png(paste(result_dir, "raw_noblacklist_overlap_rates.png", sep=""))
    #par(mfrow=c(length(unique(dbObj$samples$Condition)), 1), mar = c(5, 5, 4, 4))
    par(mar = c(5, 5, 4, 2)) 
    for (c in unique(dbObj$samples$Condition)) {
      cat('\n', c, '\n')
      olap.rate <- dba.overlap(dbObj, dba.mask(dbObj.noblacklist, attribute=DBA_CONDITION, value=c, combine='or'), mode=DBA_OLAP_RATE)
      #cat(olap.rate)
      #cat('\n')
      plot(olap.rate, type='l', xlim=c(0.75, length(olap.rate)), ylab='# overlapping peaks', xlab='# peaksets (replicates and peak callers)', col=conditions_colour_code[[c]])
      text(x=1:length(olap.rate), y=olap.rate, olap.rate, cex=0.75)
      axis(side=1, at=1:length(olap.rate))
      title(main=c("Overlap Rate - Blacklist Regions Removed", c))
    }

    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)

# Building consensus peaksets
tryCatch(
  {
    
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages

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
    cat("\n\n---------- Consensus of peak callers ----------\n\n")
    print(dbObj.caller_consensus)
    cat("\n\n---------- Final consensus from replicates ----------\n\n")
    print(dbObj.consensus)
  
    # # Add fragment sizes to new objects (helps to later run affinity analysis)
    # for (c in unique(dbObj$samples$Condition)){
    #   for (i in 1:length(dbObj.caller_consensus$mask[[c]])){
    #     if (dbObj.caller_consensus$mask[[c]][[i]]){
    #       dbObj.caller_consensus$config$fragmentSize[i] <- frag_sizes[[paste0(c, ':', i)]]
    #     }
    #   }
    # }
    dbObj.caller_consensus$config$fragmentSize <- unlist(unname(frag_sizes))
    
    for (c in unique(dbObj$samples$Condition)){
      for (i in 1:length(dbObj.consensus$mask[[c]])){
        if (dbObj.consensus$mask[[c]][[i]]){
          dbObj.consensus$config$fragmentSize[i] <- floor(mean(as.integer(frag_sizes[str_detect(names(frag_sizes), pattern=c)])))
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
  
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)


# Comparing differential/similar peaks
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    cat("\n\n---------- Analyzing Peaks ----------\n")
    
    pdf(paste(result_dir, 'binding-sites.pdf', sep=""),
        width  = figs$width,
        height = figs$height,
        pointsize = figs$pointsize)
    
    # Consensus peaks from all conditions (all relevant peaks)
    consensus_peaks <- dba.peakset(dbObj.consensus, bRetrieve=TRUE)
    
    # Create subdirectories for each condition
    result_dirs <- list()
    for (p in unique(dba.show(dbObj.consensus)$Condition)){
      result_dirs[[p]] <- paste(result_dir, p, "/", sep='')
      if (!file.exists(result_dirs[[p]])) {
        dir.create(result_dirs[[p]])
      }
    }
    result_dirs[['Shared']] <- paste(result_dir, "Shared/", sep='')
    if (!file.exists(result_dirs[["Shared"]])) {
      dir.create(result_dirs[["Shared"]])
    }
    
    # Output consensus peaksets
    raw_peaks <- list()
    for (c in unique(dba.show(dbObj.consensus)$Condition)){
      p <- dba.peakset(dbObj.consensus, dbObj.consensus$masks[[c]], bRetrieve=TRUE)
      write.table(as.data.frame(p)[c('seqnames', 'start', 'end')], file=paste(result_dirs[[c]], c, '_consensus.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
      raw_peaks[[c]] <- p
    }
    
    # Ensuring more than 1 condition provided in samplesheet
    if (length(dbObj.consensus$masks$All) > 1){
      
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
      #write.table(as.data.frame(shared_peaks[["Shared"]]), file=paste(result_dirs[[c]], '../Shared_consensus.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
      
      # If 3 conditions in samplesheet, get combinations of shared peaks
      if (length(unique(dbObj$samples$Condition)) == 3){
        for (c in unique(dba.show(dbObj.consensus)$Condition)){
          i <- which(unique(dba.show(dbObj.consensus)$Condition) == c) + 3
          pair <- unique(dba.show(dbObj.consensus)$Condition)[which(unique(dba.show(dbObj.consensus)$Condition) != c)]
          unique_peaks[[paste(pair, collapse='_and_')]] <- differential_peaks[[which(unique(dba.show(dbObj.consensus)$Condition) == c) + 3]]
          conditions_colour_code[[paste(pair, collapse='_and_')]] <- colours[i]
          cat("\n", paste(pair, collapse='_and_'), "have", length(unique_peaks[[paste(pair, collapse='_and_')]]), "shared peaks.\n")
        }
      }
      
      # Output unique and shared peaksets
      peaks <- c(unique_peaks, shared_peaks)
      
      # Output unique and shared peaks to bed files
      for (p in names(peaks)){
        if (p != "Shared"){
          write.table(as.data.frame(peaks[[p]])[c('seqnames', 'start', 'end')], file=paste(result_dirs[[p]], p, '_unique.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
        } else{
          write.table(as.data.frame(peaks[[p]])[c('seqnames', 'start', 'end')], file=paste(result_dirs[[p]], p, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
        }
      }
      
    } else {
      cat("\n\nOnly 1 condition, cannot perform a differential analysis...\n")
      for (c in unique(dba.show(dbObj.consensus)$Condition)){
        unique_peaks[[c]] <- dbObj.consensus$peaks[[1]]
        cat("\n", c, "has", dim(unique_peaks[[c]])[[1]], "unique peaks.\n")
      }
      differential_peaks <- list()
      shared_peaks <- list()
      sink()                  # stop normal output
      close(con)
      quit()
    }
    
    # # Plot peaks over genome
    # tryCatch(
    #   {
    #     plt <- covplot(c(unique_peaks[1:length(unique(dbObj$samples$Condition))], shared_peaks), title="Peaks over Genome") + #chrs = "chr1" paste0("chr", c(1:22))
    #       scale_color_manual(values=rev(c(unlist(unname(conditions_colour_code[1:length(unique(dbObj$samples$Condition))])), 'grey'))) + 
    #       scale_fill_manual(values=rev(c(unlist(unname(conditions_colour_code[1:length(unique(dbObj$samples$Condition))])), 'grey')))
    #     #invisible(capture.output(ggsave(filename=paste(result_dir, 'genome_peaks.png', sep=''), plot=plt, dpi=320)))
    #     plt <- plt + facet_grid(chr ~ .id)
    #     print(plt)
    #     #invisible(capture.output(ggsave(filename=paste(result_dir, 'genome_peaks_split.png', sep=''), plot=plt, dpi=320)))
    #     #rm(plt)
    #   },error = function(e)
    #   {
    #     message(e)
    #   }
    # )
    
    dba.plotHeatmap(dbObj.consensus, margin=15, cexRow = 0.8, cexCol = 0.8)
    #title(main=c("Correlation Heatmap - Blacklist Regions Removed"))
    mtext(c("Correlation Heatmap - Consensus Peaksets"), side = 1, line = 2)
    
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
    
    
    # Venn Plots
    if (length(unique(dba.show(dbObj.consensus)$Condition)) == 2){
      #grid.newpage()
      g = draw.pairwise.venn(area1=length(unique_peaks[[names(unique_peaks)[1]]])+length(shared_peaks[[names(shared_peaks)[1]]]), 
                             area2=length(unique_peaks[[names(unique_peaks)[2]]])+length(shared_peaks[[names(shared_peaks)[1]]]),
                             cross.area=length(shared_peaks[[names(shared_peaks)[1]]]),
                             category=names(unique_peaks),
                             fill=unname(unlist(conditions_colour_code))[1:2],
                             col=NA,
                             #cat.pos=c(0,0),
                             cat.dist = c(0,0),
                             ind=FALSE,
                             # text inside circles
                             cex = figs$pointsize/5,          # increase numbers inside areas
                             # category labels
                             cat.cex = figs$pointsize/5      # increase set names
                             #fontface = "bold" # make them bold if desired
                             )
      plt <- grid.arrange(gTree(children=g), 
                          top = textGrob("Binding Site Overlaps", gp = gpar(cex = figs$pointsize/5)),    # bigger title
                          bottom = textGrob(gsub('/', '', opt$result_dir), gp = gpar(cex = figs$pointsize/5))  # bigger bottom text
                          )
      #plot(plt)
      #invisible(capture.output(ggsave(filename=paste(result_dir, 'consensus_peaks_venn.png', sep=''), plot=plt)))
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
      plt <- plot(euler(e), main=gsub('/', '', result_dir), quantities=TRUE, fills=unname(unlist(conditions_colour_code)))
      plot(plt)
      #invisible(capture.output(ggsave(filename=paste(result_dir, 'consensus_peaks_venn.png', sep=''), plot=plt)))
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
      plt <- plot(euler(e), main=gsub('/', '', result_dir), quantities=TRUE, fills=unname(unlist(conditions_colour_code)))
      plot(plt)
      #invisible(capture.output(ggsave(filename=paste(result_dir, 'consensus_peaks_venn.png', sep=''), plot=plt)))
      #invisible(capture.output(dev.off()))
      rm(e)
    }
    
    tryCatch(
      {
        dba.plotPCA(dbObj, masks=!dbObj.total$masks$Consensus, attributes=DBA_CONDITION, label=DBA_ID, vColors=(colours_discrete),
                    labelSize  = 0.8,   # shrink or grow point labels
                    dotSize    = 1.2    # shrink or grow points
        )
        #print(plt)
        #plt$main <- "PCA"
        #invisible(capture.output(ggsave(filename=paste(result_dir, 'pca_condition.png', sep=''), plot=grid.arrange(plt))))
        
        if(length(unique(dbObj$samples$Factor)) > 1){
          dba.plotPCA(dbObj, masks=!dbObj.total$masks$Consensus, attributes=DBA_FACTOR, label=DBA_ID,
                      labelSize  = 0.8,   # shrink or grow point labels
                      dotSize    = 1.2    # shrink or grow points
          )
          #print(plt)
          #invisible(capture.output( ggsave(filename=paste(result_dir, 'pca_factor.png', sep=''), plot=grid.arrange(plt)) ))
          #invisible(capture.output(dev.off()))
        }
      },error = function(e)
      {
        message(e)
      }
    )
    
    # # Plot peaks related to TSS sites
    # tryCatch(
    #   {
    #     tagMatrices <- list()
    #     for (p in names(peaks)){
    #       tagMatrix <- getTagMatrix(peaks[[p]], windows=promoters)
    #       if (length(tagMatrix) == 0){
    #         cat("No peaks at promoter sites for", p, "\n")
    #         rm(tagMatrix)
    #         #break
    #       }else{
    #         tagMatrices[[p]] <- tagMatrix
    #         cat(dim(tagMatrix)[[1]], "peaks at promoter sites for", p, "\n")
    #         plt <- tagHeatmap(tagMatrix, 
    #                           xlab="bp at TSS", 
    #                           ylab="Peaks", 
    #                           title=paste(dim(tagMatrix)[[1]],'Peaks at Promoters', p, sep=" - "),
    #                           palette=if_else(conditions_colour_code[[p]] == "#00BFC4", 'Greens', 'Reds'), 
    #         )
    #         #invisible(capture.output(ggsave(paste(result_dirs[[p]], 'TSS_heatmap_', p, '_peaks.png', sep=''), plot=plt, dpi=320)))
    #         print(plt)
    #         rm(tagMatrix)
    #         invisible(capture.output(gc()))
    #       }
    #     }
    #   },error = function(e)
    #   {
    #     message(e)
    #   }
    # )
    # invisible(capture.output(gc()))
    
    # # Plot TSS profile of peaks
    # tryCatch(
    #   {
    #     if (length(unique(dba.show(dbObj.consensus)$Condition)) == 3){
    #       plt <- plotAvgProf(tagMatrices[1:3], xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
    #         scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[1:3]) +
    #         scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[1:3])
    #       #invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_unique-peaks.png', sep=''), plot=plt, dpi=320)))
    #       invisible(capture.output(gc()))
    #       print(plt)
    #       plt <- plotAvgProf(tagMatrices[4:6], xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
    #         scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[4:6]) +
    #         scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[4:6])
    #       #invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_pairs-peaks.png', sep=''), plot=plt, dpi=320)))
    #       invisible(capture.output(gc()))
    #       print(plt)
    #     }
    #     else{
    #       plt <- plotAvgProf(tagMatrices, xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
    #         scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))) +
    #         scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)])))
    #       #invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_peaks.png', sep=''), plot=plt, dpi=320)))
    #       invisible(capture.output(gc()))
    #       print(plt)
    #     }
    #     plt <- plotAvgProf(tagMatrices[['Shared']], xlim=c(-3000, 3000), conf=0.95, resample=1000, ncpus = parallel::detectCores()/2) +
    #       scale_color_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[length(conditions_colour_code)]) +
    #       scale_fill_manual(values=unname(unlist(conditions_colour_code[names(tagMatrices)]))[length(conditions_colour_code)])
    #     #invisible(capture.output(ggsave(paste(result_dir, 'TSS_profile_shared-peaks.png', sep=''), plot=plt, dpi=320)))
    #     invisible(capture.output(gc()))
    #     print(plt)
    #     rm(tagMatrices)
    #     rm(plt)
    #   },error = function(e)
    #   {
    #     message(e)
    #   }
    # )
    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)


# Annotating consensus peaks
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    cat("\n\n---------- Annotating Consensus Peaks ----------\n")
    
    pdf(paste(result_dir, 'annotated-raw-consensus-sites.pdf', sep=""),
        width  = figs$width,
        height = figs$height,
        pointsize = figs$pointsize)  #8 or 9–10)
    
    # ========= Get Annotations =========
    peakAnnoList <- list()
      
    # Annotate consensus peaks
    for (p in names(raw_peaks)){
      # Set plot colors
      colour <- conditions_colour_code[[p]]
      heat_colour <- conditions_colour_code[[p]]
      
      cat("\nAnnotating", p, 'Consensus Peaks\n')
      anno <- annotatePeak(raw_peaks[[p]], 
                           TxDb=anno_ref$txdb,
                           annoDb=anno_ref$annoDb,
                           level=opt$annotation_level,
                           tssRegion=c(-3000, 3000)
      )
      peakAnnoList[[p]] <- anno
      
      cat("\n",length(anno@anno), "annotated out of", length(raw_peaks[[p]]), p, "peaks\n")
      plt <- make_anno_piebar(as.data.frame(anno@anno), type='pie', title=paste0(p, " (Consensus)\n", "Distribution of Sites"), specific=TRUE, colours=colours_discrete, 
                              title_size=figs$width*2.8,
                              text_size=figs$pointsize/3,
                              legend_key_size=figs$width/20,
                              legend_title_size=figs$width*1.8,
                              legend_text_size=figs$width*1.5)
      print(plt)
      plt <- make_anno_piebar(as.data.frame(anno@anno), type='bar', title=paste0(p, " (Consensus)\n", "Distribution of Sites"), specific=TRUE, colours=colours_discrete,
                              title_size=figs$width*2,
                              text_size=figs$pointsize/2.8,
                              axis_title_size=figs$width*1.6,
                              axis_x_size=figs$width*1.3,
                              axis_y_size=figs$width*1.3
                              )
      print(plt)
      write.table(anno@anno, file=paste(result_dirs[[p]], p, '_consensus_annotated.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
      
      
      # Mapper for EntrezID to gene SYMBOL
      mapper <- as.data.frame(anno@anno)
      mapper <- mapper[c('geneId', 'SYMBOL')]
      mapper <- mapper[!duplicated(mapper), ]
      
      tryCatch(
        {
          genes <- list()
          genes[[p]] <- anno@anno$geneId
          names(genes) = sub("_", "\n", names(genes))
          
          # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
          compKEGG <- compareCluster(geneCluster=genes,
                                     fun="enrichKEGG",
                                     pvalueCutoff=opt$fdr,
                                     pAdjustMethod="BH",
                                     minGSSize = opt$minGSSize,
                                     maxGSSize = opt$maxGSSize,
                                     organism=anno_ref$keggOrg,
                                     universe = bg_entrez
          ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          if (class(compKEGG) == 'compareClusterResult'){
            if (dim(compKEGG@compareClusterResult)[1] > 0) {
              cat('\n', dim(compKEGG@compareClusterResult)[1], 'KEGG results\n')
              # Map EntrezIDs to gene SYMBOL
              compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
              myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
              for (i in 1:length(myEntrez)){
                compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
              }
              # Write annotations to csv
              df_kegg <- as.data.frame(compKEGG@compareClusterResult)
              df_kegg <- df_kegg[order(df_kegg$p.adjust, -xtfrm(df_kegg$GeneRatio), -xtfrm(df_kegg$BgRatio), df_kegg$Description), ]
              write.table(df_kegg, file=paste(result_dirs[[p]], p, '_consensus_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_anno_dotplot(df_kegg, 
                                       title=paste('KEGG - ', p,  ' (Consensus)', sep=""), 
                                       ylabel="KEGG Category", 
                                       colour=colours_continuous, 
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
            }
            
          } else{
            cat('\nNo KEGG results.\n')
          }
          
        },error = function(e)
        {
          message(e)
          gc()
        }
      )
      invisible(capture.output(gc()))
      
      
      for (ont in c('ALL', 'CC', 'MF', 'BP')){
        tryCatch(
          {
            genes <- list()
            genes[[p]] <- anno@anno$SYMBOL
            names(genes) = sub("_", "\n", names(genes))
            #cat("\nGetting ", p, ' GO ', ont, '\n')
            
            compGO <- compareCluster(geneCluster=genes,
                                     keyType='SYMBOL',
                                     OrgDb=anno_ref$annoDb,
                                     fun="enrichGO",
                                     ont=ont,
                                     pvalueCutoff=opt$fdr,
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
              df_go <- as.data.frame(compGO@compareClusterResult)
              df_go <- df_go[order(df_go$p.adjust, -xtfrm(df_go$GeneRatio), -xtfrm(df_go$BgRatio), df_go$Description), ]
              write.table(df_go, file=paste(result_dirs[[p]], p, '_consensus_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_anno_dotplot(df_go, 
                                       title=paste("GO (", ont, ") - ", p, " (Consensus)", sep=""), 
                                       ylabel="GO Term", 
                                       colour=colours_continuous, 
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
              gc()
            } else{
              cat("\nNo GO", ont, "annotation results\n")
              gc()
            }
          },error = function(e)
          {
            message(e)
            gc()
          }
        )
        invisible(capture.output(gc()))
      }
      
      
      genes_entrez <- list()
      genes_entrez[[p]] <- anno@anno$geneId
      names(genes_entrez) = sub("_", "\n", names(genes_entrez))
      for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
        #cat("\nGetting ", p, ' DAVID ', annotation_type, ' annotations...\n')  
        # DAVID Annotation
        tryCatch(
          {
            #cat("\nDAVID - ", annotation_type, "\n")
            entrez <- unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))])
            
            py_require(c("pandas"))
            py_require(c("suds"))
            py_run_string("import pandas as pd")
            py_run_string("import sys")
            py_run_string("from suds.client import Client")
            
            # create a service client using the wsdl.
            py_run_string("client = Client('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl')")
            py_run_string("client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')")
            
            #authenticate user email
            py_run_string("client.service.authenticate(r.opt['david_user'])")
            
            # Read input gene list file, convert ids to a comma-delimited string and upload the list to DAVID
            py_run_string("inputIds = ','.join(r.entrez)")
            py_run_string("client.service.addList(inputIds, 'ENTREZ_GENE_ID', r.p, 0)")
            
            # setCategories
            py_run_string("categorySting = str(client.service.setCategories(r.annotation_type))")
            
            #getChartReport
            py_run_string("thd = r.opt['fdr']")
            py_run_string("ct = r.opt['minGSSize']")
            py_run_string("chartReport = client.service.getChartReport(thd,ct)")
            py_run_string("chartRow = len(chartReport)")
            py_run_string("print ('Total chart records:',chartRow)")
          },error = function(e)
          {
            cat("\nUnable to run DAVID\n")
            message(e)
            next
            gc()
          }
        )
        tryCatch(
          {
            if (py$chartRow > 0){
              if (annotation_type == "KEGG_PATHWAY"){
                splitter <- ":"
              }else{
                splitter <- "~"
              }
              # parse chartReport
              records <- data.frame(
                ID            = character(),
                Category      = character(),
                Description   = character(),
                GeneRatio     = character(),
                BgRatio       = character(),
                pvalue        = numeric(),
                p.adjust      = numeric(),
                FDR           = numeric(),
                geneID        = character(),
                Count         = integer(),
                foldEnrichment= numeric(),
                id            = character(),
                stringsAsFactors = FALSE
              )
              
              for (simpleChartRecord in py$chartReport) {
                df_record <- data.frame(
                  ID            = strsplit(simpleChartRecord$termName, splitter)[[1]][1],
                  Category      = simpleChartRecord$categoryName,
                  Description   = strsplit(simpleChartRecord$termName, splitter)[[1]][2],
                  GeneRatio     = paste0(simpleChartRecord$listHits, "/", simpleChartRecord$listTotals),
                  BgRatio       = paste0(simpleChartRecord$popHits, "/", simpleChartRecord$popTotals),
                  pvalue        = simpleChartRecord$ease,
                  p.adjust      = simpleChartRecord$benjamini,
                  FDR           = simpleChartRecord$afdr,
                  geneID        = gsub(", ", "/", simpleChartRecord$geneIds),
                  Count         = simpleChartRecord$listHits,
                  foldEnrichment= simpleChartRecord$foldEnrichment,
                  id            = simpleChartRecord$id,
                  stringsAsFactors = FALSE
                )
                records <- rbind(records, df_record)
              }
              # py_run_string("records = pd.DataFrame()")
              # py_run_string("for simpleChartRecord in chartReport:
              #                 df = pd.DataFrame(index=[records.shape[0]], data={
              #                     'ID' : simpleChartRecord['termName'].split(r.splitter)[0],
              #                     'Category' : simpleChartRecord['categoryName'],
              #                     'Description' : simpleChartRecord['termName'].split(r.splitter)[1],
              #                     'GeneRatio': str(simpleChartRecord['listHits']) + '/' + str(simpleChartRecord['listTotals']), 
              #                     'BgRatio': str(simpleChartRecord['popHits']) + '/' + str(simpleChartRecord['popTotals']), 
              #                     'pvalue' : simpleChartRecord['ease'],
              #                     'p.adjust' : simpleChartRecord['benjamini'],
              #                     'FDR' : simpleChartRecord['afdr'],
              #                     'geneID' : simpleChartRecord['geneIds'].replace(', ', '/'),
              #                     'Count' : simpleChartRecord['listHits'],
              #                     'foldEnrichment' : simpleChartRecord['foldEnrichment'],
              #                     'id' : simpleChartRecord['id']
              #                     })
              #                 records = pd.concat([records, df])")
              # #py_run_string("records.reset_index(drop=True, inplace=True)")
              #compD <- py$records
              compD <- records
              colnames(compD)[colnames(compD) == "foldEnrichment"] <- "FoldEnrichment"
              qobj <- tryCatch(qvalue(p=as.numeric(compD$pvalue), lambda=opt$pvalue, pi0.method="bootstrap"),
                               error=function(e) NULL)
              if (class(qobj) == "qvalue") {
                qvalues <- qobj$qvalues
              } else {
                qvalues <- NA
              }
              compD$qvalue <- qvalues
              
              compDAVID <- new("enrichResult",
                               result         = compD,
                               pvalueCutoff   = opt$fdr,
                               pAdjustMethod  = "BH",
                               organism       = opt$assembly,
                               ontology       = annotation_type,
                               gene           = entrez,
                               keytype        = "ENTREZ_GENE_ID")
              #rm(compD)
              
              if (!class(compDAVID) == 'enrichResult'){
                cat("\nNo DAVID", annotation_type, "results.\n")
                next
              }else{
                cat('\n', dim(compDAVID@result)[1], 'DAVID', annotation_type,  'results\n')
              }
              # old deprecated function
              #compDAVID <- enrichDAVID(
              #  unname(genes_entrez[[n]][!is.na(unname(genes_entrez[[n]]))]),
              #  idType = "ENTREZ_GENE_ID",
              #  minGSSize = opt$minGSSize,
              #  maxGSSize = opt$maxGSSize,
              #  annotation = annotation_type,
              #  pvalueCutoff = opt$pvalue,
              #  pAdjustMethod = "BH",
              #  #species = NA,
              #  david.user=opt$david_user
              #)
              
              if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
                # Map EntrezIDs to gene SYMBOL
                compDAVID@result$SYMBOL <- compDAVID@result$geneID
                myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
                for (i in 1:length(myEntrez)){
                  compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                }
                # Write annotations to csv
                df_david <- as.data.frame(compDAVID@result)
                df_david <- df_david[order(df_david$p.adjust, df_david$FoldEnrichment, -xtfrm(df_david$GeneRatio), -xtfrm(df_david$BgRatio), df_david$Description), ]
                write.table(df_david, file=paste(result_dirs[[p]], p, '_consensus_annotated_DAVID_', annotation_type, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_anno_dotplot(df_david, 
                                         title=paste('DAVID - ', p, " (Consensus)", sep=""), 
                                         ylabel=paste(annotation_type,"Category", sep=' '), 
                                         colour=colours_continuous, 
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
                #invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], 'DAVID_annotation_', annotation_type, '_', p, '_dotplot.png', sep=''), plot=plt, dpi=320)))
                print(plt)
                #remove(compDAVID)
                gc()
              } else{
                next
                cat("\nNo DAVID", annotation_type, "annotation results\n")
                gc()
              }
            }
          },error = function(e)
          {
            message(e)
            next
            gc()
          }
        )
      }
      
      
      
    }
    
    # Plot peak distributions
    tryCatch(
      {
        plt <- plotAnnoBar(peakAnnoList, title = 'Feature Distribution of Consensus Peaks')
        plt + scale_fill_manual(values = colours_discrete)
        print(plt)
        #invisible(capture.output(ggsave(filename=paste(result_dir, 'peaks_annotation_distribution_bar.png', sep=''), plot=plt, dpi=320)))
      },error = function(e)
      {
        message(e)
      }
    )
    tryCatch(
      {
        plt <- plotDistToTSS(peakAnnoList, title = 'Feature Distribution of Consensus Peaks Relative to TSS')
        plt + scale_fill_manual(values = colours_discrete)
        print(plt)
        #invisible(capture.output(ggsave(filename=paste(result_dir, 'peaks_annotation_TSS_distribution.png', sep=''), plot=plt, dpi=320)))
      },error = function(e)
      {
        message(e)
      }
    )
    invisible(capture.output(gc()))
    
    # # Annotate shared peaks
    # if (length(dbObj.consensus$masks$All) > 1){
    #   # Set plot colors
    #   colour <- "#56B1F7"
    #   heat_colour <- "BuPu"
    #   anno <- annotatePeak(shared_peaks[["Shared"]], 
    #                        TxDb=anno_ref$txdb,
    #                        annoDb=anno_ref$annoDb,
    #                        level=opt$annotation_level,
    #                        tssRegion=c(-3000, 3000)
    #   )
    #   cat("\n", length(anno@anno), "annotated out of", length(shared_peaks[['Shared']]), "shared peaks\n\n")
    #   plt <- make_anno_piebar(as.data.frame(anno@anno), type='pie', title=paste0("Shared\n", "Distribution of Sites"), specific=TRUE, colours=paletteer_d("khroma::muted"))
    #   print(plt)
    #   plt <-make_anno_piebar(as.data.frame(anno@anno), type='bar', title=paste0("Shared\n", "Distribution of Sites"), specific=TRUE, colours=paletteer_d("khroma::muted"))
    #   print(plt)
    #   write.table(anno@anno, file=paste(result_dirs[[c]], '../Shared_annotated.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
    
    if (!length(dbObj.consensus$masks$All) > 1) {
      cat("\nnOnly 1 condition, cannot perform a differential analysis...\n")
      for (c in unique(dba.show(dbObj.consensus)$Condition)){
        unique_peaks[[c]] <- dbObj.consensus$peaks[[1]]
        cat("\n", c, "has", dim(unique_peaks[[c]])[[1]], "unique peaks.\n")
      }
      differential_peaks <- list()
      shared_peaks <- list()
      sink()                  # stop normal output
      close(con)
      quit()
    }
    
    # Upset plot for annotated genes
    upsetlist <- list()
    for (p in names(peakAnnoList)) {
      upsetlist[[p]] <- peakAnnoList[[p]]@anno$SYMBOL
    }
    upset_colors <- list()
    for (n in names(peakAnnoList)){
      upset_colors[[ conditions_colour_code[[n]] ]] <- length(unique(peakAnnoList[[n]]@anno$SYMBOL))
    }
    upset_colors <- sort(unlist(upset_colors), decreasing=TRUE)
    
    plt <- upset(fromList(upsetlist), 
          order.by = "freq", 
          nsets = length(names(peakAnnoList)),
          sets.bar.color = names(upset_colors),
          empty.intersections = "on",
          set_size.show = TRUE,
          set_size.angles = 0,
          set_size.scale_max = dim(fromList(upsetlist))[[1]]*1.25,
          sets.x.label = "Gene Set Size",
          mainbar.y.label = "Intersection Size of Gene Sets",
          mb.ratio = c(0.65, 0.35),
          #main.bar.color = unlist(conditions_colour_code),
          text.scale = c(figs$pointsize/6, figs$pointsize/6, figs$pointsize/7, figs$pointsize/7, figs$pointsize/6, figs$pointsize/6)
          # order: intersection size title, set size title,
          # intersection tick labels, set tick labels,
          # intersection bar labels, set size labels
    )
    print(plt)
    grid.text("Annotated Genes in Consensus Peaksets",x = 0.65, y=0.95, gp=gpar(fontsize=figs$pointsize*1.2))
    

  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)


# Annotating differential and shared peaks
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    cat("\n\n---------- Annotating Differential and Shared Peaks ----------\n")
    
    pdf(paste(result_dir, 'annotated-differential-and-shared-sites.pdf', sep=""),
        width  = figs$width,
        height = figs$height,
        pointsize = figs$pointsize)  #8 or 9–10)
    
    # ========= Get Annotations =========
    peakAnnoList <- list()
    
    # Annotate consensus peaks
    for (p in names(peaks)){
      # Set plot colors
      colour <- conditions_colour_code[[p]]
      heat_colour <- conditions_colour_code[[p]]
      
      cat("\nAnnotating", p, 'Unique Peaks\n')
      anno <- annotatePeak(peaks[[p]], 
                           TxDb=anno_ref$txdb,
                           annoDb=anno_ref$annoDb,
                           level=opt$annotation_level,
                           tssRegion=c(-3000, 3000)
      )
      peakAnnoList[[p]] <- anno
      
      cat("\n",length(anno@anno), "annotated out of", length(peaks[[p]]), p, "peaks\n")
      plt <- make_anno_piebar(as.data.frame(anno@anno), type='pie', title=paste0(p, " - Unique\n", "Distribution of Sites"), specific=TRUE, colours=paletteer_d("khroma::muted"), 
                              title_size=figs$width*2.8,
                              text_size=figs$pointsize/3,
                              legend_key_size=figs$width/20,
                              legend_title_size=figs$width*1.8,
                              legend_text_size=figs$width*1.5)
      print(plt)
      plt <- make_anno_piebar(as.data.frame(anno@anno), type='bar', title=paste0(p, " - Unique\n", "Distribution of Sites"), specific=TRUE, colours=paletteer_d("khroma::muted"),
                              title_size=figs$width*2,
                              text_size=figs$pointsize/2.8,
                              axis_title_size=figs$width*1.6,
                              axis_x_size=figs$width*1.3,
                              axis_y_size=figs$width*1.3
      )
      print(plt)
      write.table(anno@anno, file=paste(result_dirs[[p]], p, '_unique_annotated.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
      
      
      # Mapper for EntrezID to gene SYMBOL
      mapper <- as.data.frame(anno@anno)
      mapper <- mapper[c('geneId', 'SYMBOL')]
      mapper <- mapper[!duplicated(mapper), ]
      
      tryCatch(
        {
          genes <- list()
          genes[[p]] <- anno@anno$geneId
          names(genes) = sub("_", "\n", names(genes))
          
          # fun is "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" 
          compKEGG <- compareCluster(geneCluster=genes,
                                     fun="enrichKEGG",
                                     pvalueCutoff=opt$fdr,
                                     pAdjustMethod="BH",
                                     minGSSize = opt$minGSSize,
                                     maxGSSize = opt$maxGSSize,
                                     organism=anno_ref$keggOrg,
                                     universe = bg_entrez
          ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
          if (class(compKEGG) == 'compareClusterResult'){
            
            if (dim(compKEGG@compareClusterResult)[1] > 0){
              cat('\n', dim(compKEGG@compareClusterResult)[1], 'KEGG results\n')
              # Map EntrezIDs to gene SYMBOL
              compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
              myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
              for (i in 1:length(myEntrez)){
                compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
              }
              # Write annotations to csv
              df_kegg <- as.data.frame(compKEGG@compareClusterResult)
              df_kegg <- df_kegg[order(df_kegg$p.adjust, -xtfrm(df_kegg$GeneRatio), -xtfrm(df_kegg$BgRatio), df_kegg$Description), ]
              write.table(df_kegg, file=paste(result_dirs[[p]], p, '_unique_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_anno_dotplot(df_kegg, 
                                       title=paste('KEGG - ', p,  ' Unique', sep=""), 
                                       ylabel="KEGG Category", 
                                       #colour=colour, 
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
              
            }
            
          } else{
              cat("\nNo KEGG results.\n")
              #next
          }
          
        },error = function(e)
        {
          message(e)
          gc()
        }
      )
      invisible(capture.output(gc()))
      
      
      for (ont in c('ALL', 'CC', 'MF', 'BP')){
        tryCatch(
          {
            genes <- list()
            genes[[p]] <- anno@anno$SYMBOL
            names(genes) = sub("_", "\n", names(genes))
            #cat("\nGetting ", p, ' GO ', ont, '\n')
            
            compGO <- compareCluster(geneCluster=genes,
                                     keyType='SYMBOL',
                                     OrgDb=anno_ref$annoDb,
                                     fun="enrichGO",
                                     ont=ont,
                                     pvalueCutoff=opt$fdr,
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
              df_go <- as.data.frame(compGO@compareClusterResult)
              df_go <- df_go[order(df_go$p.adjust, -xtfrm(df_go$GeneRatio), -xtfrm(df_go$BgRatio), df_go$Description), ]
              write.table(df_go, file=paste(result_dirs[[p]], p, '_unique_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
              
              plt <- make_anno_dotplot(df_go, 
                                       title=paste("GO (", ont, ") - ", p, " Consensus", sep=""), 
                                       ylabel="GO Term", 
                                       #colour=colour, 
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
              gc()
            } else{
              cat("\nNo GO", ont, " results\n")
              gc()
            }
          },error = function(e)
          {
            message(e)
            gc()
          }
        )
        invisible(capture.output(gc()))
      }
      
      
      genes_entrez <- list()
      genes_entrez[[p]] <- anno@anno$geneId
      names(genes_entrez) = sub("_", "\n", names(genes_entrez))
      for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
        #cat("\nGetting ", p, ' DAVID ', annotation_type, ' annotations...\n')  
        # DAVID Annotation
        tryCatch(
          {
            #cat("\nDAVID - ", annotation_type, "\n")
            entrez <- unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))])
            
            py_require(c("pandas"))
            py_require(c("suds"))
            py_run_string("import pandas as pd")
            py_run_string("import sys")
            py_run_string("from suds.client import Client")
            
            # create a service client using the wsdl.
            py_run_string("client = Client('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl')")
            py_run_string("client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')")
            
            #authenticate user email
            py_run_string("client.service.authenticate(r.opt['david_user'])")
            
            # Read input gene list file, convert ids to a comma-delimited string and upload the list to DAVID
            py_run_string("inputIds = ','.join(r.entrez)")
            py_run_string("client.service.addList(inputIds, 'ENTREZ_GENE_ID', r.p, 0)")
            
            # setCategories
            py_run_string("categorySting = str(client.service.setCategories(r.annotation_type))")
            
            #getChartReport
            py_run_string("thd = r.opt['fdr']")
            py_run_string("ct = r.opt['minGSSize']")
            py_run_string("chartReport = client.service.getChartReport(thd,ct)")
            py_run_string("chartRow = len(chartReport)")
            py_run_string("print ('Total chart records:',chartRow)")
            
          },error = function(e)
          {
            cat("\nUnable to run DAVID\n")
            message(e)
            next
            gc()
          }
        )
        tryCatch(
          {   
            if (py$chartRow > 0){
              if (annotation_type == "KEGG_PATHWAY"){
                splitter <- ":"
              }else{
                splitter <- "~"
              }
              # parse chartReport
              records <- data.frame(
                ID            = character(),
                Category      = character(),
                Description   = character(),
                GeneRatio     = character(),
                BgRatio       = character(),
                pvalue        = numeric(),
                p.adjust      = numeric(),
                FDR           = numeric(),
                geneID        = character(),
                Count         = integer(),
                foldEnrichment= numeric(),
                id            = character(),
                stringsAsFactors = FALSE
              )
              
              for (simpleChartRecord in py$chartReport) {
                df_record <- data.frame(
                  ID            = strsplit(simpleChartRecord$termName, splitter)[[1]][1],
                  Category      = simpleChartRecord$categoryName,
                  Description   = strsplit(simpleChartRecord$termName, splitter)[[1]][2],
                  GeneRatio     = paste0(simpleChartRecord$listHits, "/", simpleChartRecord$listTotals),
                  BgRatio       = paste0(simpleChartRecord$popHits, "/", simpleChartRecord$popTotals),
                  pvalue        = simpleChartRecord$ease,
                  p.adjust      = simpleChartRecord$benjamini,
                  FDR           = simpleChartRecord$afdr,
                  geneID        = gsub(", ", "/", simpleChartRecord$geneIds),
                  Count         = simpleChartRecord$listHits,
                  foldEnrichment= simpleChartRecord$foldEnrichment,
                  id            = simpleChartRecord$id,
                  stringsAsFactors = FALSE
                )
                records <- rbind(records, df_record)
              }
              # py_run_string("records = pd.DataFrame()")
              # py_run_string("for simpleChartRecord in chartReport:
              #                 df = pd.DataFrame(index=[records.shape[0]], data={
              #                     'ID' : simpleChartRecord['termName'].split(r.splitter)[0],
              #                     'Category' : simpleChartRecord['categoryName'],
              #                     'Description' : simpleChartRecord['termName'].split(r.splitter)[1],
              #                     'GeneRatio': str(simpleChartRecord['listHits']) + '/' + str(simpleChartRecord['listTotals']), 
              #                     'BgRatio': str(simpleChartRecord['popHits']) + '/' + str(simpleChartRecord['popTotals']), 
              #                     'pvalue' : simpleChartRecord['ease'],
              #                     'p.adjust' : simpleChartRecord['benjamini'],
              #                     'FDR' : simpleChartRecord['afdr'],
              #                     'geneID' : simpleChartRecord['geneIds'].replace(', ', '/'),
              #                     'Count' : simpleChartRecord['listHits'],
              #                     'foldEnrichment' : simpleChartRecord['foldEnrichment'],
              #                     'id' : simpleChartRecord['id']
              #                     })
              #                 records = pd.concat([records, df])")
              # #py_run_string("records.reset_index(drop=True, inplace=True)")
              #compD <- py$records
              compD <- records
              colnames(compD)[colnames(compD) == "foldEnrichment"] <- "FoldEnrichment"
              qobj <- tryCatch(qvalue(p=as.numeric(compD$pvalue), lambda=opt$pvalue, pi0.method="bootstrap"),
                               error=function(e) NULL)
              if (class(qobj) == "qvalue") {
                qvalues <- qobj$qvalues
              } else {
                qvalues <- NA
              }
              compD$qvalue <- qvalues
              
              compDAVID <- new("enrichResult",
                               result         = compD,
                               pvalueCutoff   = opt$fdr,
                               pAdjustMethod  = "BH",
                               organism       = opt$assembly,
                               ontology       = annotation_type,
                               gene           = entrez,
                               keytype        = "ENTREZ_GENE_ID")
              #rm(compD)
              
              if (!class(compDAVID) == 'enrichResult'){
                cat("\nNo DAVID", annotation_type, "results.\n")
                next
              }else{
                cat('\n', dim(compDAVID@result)[1], 'DAVID', annotation_type,  'results\n')
              }
              # old deprecated function
              #compDAVID <- enrichDAVID(
              #  unname(genes_entrez[[n]][!is.na(unname(genes_entrez[[n]]))]),
              #  idType = "ENTREZ_GENE_ID",
              #  minGSSize = opt$minGSSize,
              #  maxGSSize = opt$maxGSSize,
              #  annotation = annotation_type,
              #  pvalueCutoff = opt$pvalue,
              #  pAdjustMethod = "BH",
              #  #species = NA,
              #  david.user=opt$david_user
              #)
              
              if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
                # Map EntrezIDs to gene SYMBOL
                compDAVID@result$SYMBOL <- compDAVID@result$geneID
                myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
                for (i in 1:length(myEntrez)){
                  compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                }
                # Write annotations to csv
                df_david <- as.data.frame(compDAVID@result)
                df_david <- df_david[order(df_david$p.adjust, df_david$FoldEnrichment, -xtfrm(df_david$GeneRatio), -xtfrm(df_david$BgRatio), df_david$Description), ]
                write.table(df_david, file=paste(result_dirs[[p]], p, '_unique_annotated_DAVID_', annotation_type, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_anno_dotplot(df_david, 
                                         title=paste('DAVID - ', p, " Unique", sep=""), 
                                         ylabel=paste(annotation_type,"Category", sep=' '), 
                                         #colour=colour, 
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
                #invisible(capture.output(ggsave(filename=paste(result_dirs[[p]], 'DAVID_annotation_', annotation_type, '_', p, '_dotplot.png', sep=''), plot=plt, dpi=320)))
                print(plt)
                gc()
              } else{
                cat("\nNo DAVID", annotation_type, "annotation results\n")
                next
                gc()
              }
            }
          },error = function(e)
          {
            message(e)
            next
            gc()
          }
        )
      }
      
      
      
    }
    
    # Plot peak distributions
    tryCatch(
      {
        plt <- plotAnnoBar(peakAnnoList, title = 'Feature Distribution of Unique Peaks')
        print(plt)
        #invisible(capture.output(ggsave(filename=paste(result_dir, 'peaks_annotation_distribution_bar.png', sep=''), plot=plt, dpi=320)))
      },error = function(e)
      {
        message(e)
      }
    )
    tryCatch(
      {
        plt <- plotDistToTSS(peakAnnoList, title = 'Feature Distribution of Unique Peaks Relative to TSS')
        print(plt)
        #invisible(capture.output(ggsave(filename=paste(result_dir, 'peaks_annotation_TSS_distribution.png', sep=''), plot=plt, dpi=320)))
      },error = function(e)
      {
        message(e)
      }
    )
    invisible(capture.output(gc()))
    
    # # Annotate shared peaks
    # if (length(dbObj.consensus$masks$All) > 1){
    #   # Set plot colors
    #   colour <- "#56B1F7"
    #   heat_colour <- "BuPu"
    #   anno <- annotatePeak(shared_peaks[["Shared"]], 
    #                        TxDb=anno_ref$txdb,
    #                        annoDb=anno_ref$annoDb,
    #                        level=opt$annotation_level,
    #                        tssRegion=c(-3000, 3000)
    #   )
    #   cat("\n", length(anno@anno), "annotated out of", length(shared_peaks[['Shared']]), "shared peaks\n\n")
    #   plt <- make_anno_piebar(as.data.frame(anno@anno), type='pie', title=paste0("Shared\n", "Distribution of Sites"), specific=TRUE, colours=paletteer_d("khroma::muted"))
    #   print(plt)
    #   plt <-make_anno_piebar(as.data.frame(anno@anno), type='bar', title=paste0("Shared\n", "Distribution of Sites"), specific=TRUE, colours=paletteer_d("khroma::muted"))
    #   print(plt)
    #   write.table(anno@anno, file=paste(result_dirs[[c]], '../Shared_annotated.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
    
    if (!length(dbObj.consensus$masks$All) > 1) {
      cat("\nnOnly 1 condition, cannot perform a differential analysis...\n")
      for (c in unique(dba.show(dbObj.consensus)$Condition)){
        unique_peaks[[c]] <- dbObj.consensus$peaks[[1]]
        cat("\n", c, "has", dim(unique_peaks[[c]])[[1]], "unique peaks.\n")
      }
      differential_peaks <- list()
      shared_peaks <- list()
      sink()                  # stop normal output
      close(con)
      quit()
    }
    
    # Upset plot for annotated genes
    upsetlist <- list()
    for (p in names(peakAnnoList[ names(peakAnnoList) != "Shared" ])) {
      upsetlist[[p]] <- peakAnnoList[[p]]@anno$SYMBOL
    }
    upset_colors <- list()
    for (n in names(peakAnnoList[ names(peakAnnoList) != "Shared" ])){
      upset_colors[[ conditions_colour_code[[n]] ]] <- length(unique(peakAnnoList[[n]]@anno$SYMBOL))
    }
    upset_colors <- sort(unlist(upset_colors), decreasing=TRUE)
    
    plt <- upset(fromList(upsetlist), 
                 order.by = "freq", 
                 nsets = length(names(peakAnnoList[ names(peakAnnoList) != "Shared" ])),
                 sets.bar.color = names(upset_colors),
                 empty.intersections = "on",
                 set_size.show = TRUE,
                 set_size.angles = 0,
                 set_size.scale_max = dim(fromList(upsetlist))[[1]]*1.25,
                 sets.x.label = "Gene Set Size",
                 mainbar.y.label = "Intersection Size of Gene Sets",
                 mb.ratio = c(0.65, 0.35),
                 text.scale = c(figs$pointsize/6, figs$pointsize/6, figs$pointsize/7, figs$pointsize/7, figs$pointsize/6, figs$pointsize/6)
                 # order: intersection size title, set size title,
                 # intersection tick labels, set tick labels,
                 # intersection bar labels, set size labels
    )
    print(plt)
    grid.text("Annotated Genes in Unique Peaksets",x = 0.65, y=0.95, gp=gpar(fontsize=figs$pointsize*1.2))
    
    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)

# ========= END OF OCCUPANCY ANALYSIS =========
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages

    cat("\n\n=============================== END OF OCCUPANCY ANALYSIS ===============================\n")
    if (opt$occupancy_only){
      cat("Skipping affinity analysis...\nDone!\n")
      sink()                  # stop normal output
      close(con)
      q()
    }
    # Free up memory
    rm(dbObj.final)
    rm(dbObj.total)
    rm(plt)
    rm(anno)
    rm(peakAnnoList)
    rm(compD)
    rm(compDAVID)
    rm(compKEGG)
    rm(compGO)
    rm(records)
    rm(df_record)
    rm(df_david)
    rm(df_kegg)
    rm(df_go)
    invisible(capture.output(gc()))
    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)

# ========= START OF AFFINITY ANALYSIS =========
result_dir <- paste(opt$result_dir, 'Affinity_Analysis/', sep='')
# Make new directory for analysis
if (!file.exists(result_dir)) {
  dir.create(result_dir)
}

if (!file.exists(paste(result_dir, 'DESeq2/', sep=''))) {
  dir.create(paste(result_dir, 'DESeq2/', sep=''))
}
if (!file.exists(paste(result_dir, 'edgeR/', sep=''))) {
  dir.create(paste(result_dir, 'edgeR/', sep=''))
}

output_prefix <- change_dirs(result_dir, '', '')

tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    
    cat(
      "\n\n=============================== START OF AFFINITY ANALYSIS ===============================
    The consensus peaks determined in the previous occupancy analysis are used here to focus
    on relevant peak regions only. Read counts for those peaks are analyzed for statistical 
    significance in differential binding comparisons. 
    DESeq2, edgeR, and overlapping results are produced."
    )

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
    } else{
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
    #cat("\n\n---------- Reads in Peaks Counted ----------\n")
    #print(dbObj.counted)
    
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
        message(e, "\nUsing approximate normalization...\n")
        # Approximate normalization to above without extra reading of bam files
        dbObj.norm <- dba.normalize(dbObj.counted, method=DBA_ALL_METHODS, 
                                    normalize=DBA_NORM_LIB,
                                    background=FALSE, library=DBA_LIBSIZE_FULL,
                                    spikein=FALSE, offsets=FALSE,
                                    libFun=mean, bRetrieve=FALSE)
      }
    )
    cat("\n\n---------- Reads in Peaks Counted and Normalized ----------\n\n")
    print(dbObj.norm)
    
  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)


# ========= Define contrasts between sample conditions (assumes 2 conditions) =========
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages

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
    cat("\n\n---------- Contrast Defined ----------\n\n")
    print(dbObj.contrast)

    # ========= Analyze Affinities =========
    dbObj.analyzed <- dba.analyze(dbObj.contrast, method=DBA_ALL_METHODS, bParallel=TRUE)
    cat("\n\n---------- Analyzed ----------\n\n")
    print(dbObj.analyzed)

  },
  error = function(e) {
    message("Error occurred: ", conditionMessage(e))
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)

# Generate Plots
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages
    
    pdf(paste(result_dir, 'affinity-analysis_binding-sites.pdf', sep=""),
        width  = figs$width,
        height = figs$height,
        pointsize = figs$pointsize)  # or 9–10)

    #png(paste(result_dir, 'consensus_peaks_counted_normalized_heatmap.png', sep=''))
    dba.plotHeatmap(dbObj.norm, margin=15, cexRow = 0.8, cexCol = 0.8)
    mtext(c("Correlation Heatmap - Counts Normalized"), side = 1, line = 2)

    dba.plotPCA(dbObj.norm, attributes=DBA_CONDITION, label=DBA_ID, vColors=(colours_discrete),
                    labelSize  = 0.8,   # shrink or grow point labels
                    dotSize    = 1.2    # shrink or grow points
                    )

    
    tryCatch(
      {
        dba.plotVenn(main="Differential Binding Sites Identified by Method",
                            dbObj.analyzed, 
                            contrast=1, 
                            method=DBA_ALL_METHODS
        )
      }, error=function(e){
        message("\nNo venn diagram for number of DB sites by each method\n", e)
      }
    )
    
    for (m in c(DBA_DESEQ2, DBA_EDGER)){
      tryCatch(
        {

          dba.plotVenn(dbObj.analyzed, method=m, contrast=1, 
                        bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE,
                        main=paste(dbObj.analyzed$contrasts[[1]]$name2, " Binding Sites Identified by ", m, sep='')
                      )
        }, error=function(e){
          message('\nNo venn of differential peaks from', m, '\n')
        }
      )

      tryCatch(
        {
          dba.plotVolcano(dbObj.analyzed, method=m)
        }, error=function(e){
          message('\nNo volcano plot of differential peaks from', m, '\n')
        }
      )

      tryCatch(
        {
          dba.plotMA(dbObj.analyzed, method=m)
        }, error=function(e){
          message('\nNo MA plot of differential peaks from', m, '\n')
        }
      )

      tryCatch(
        {
          dba.plotBox(dbObj.analyzed, method=m)
        }, error=function(e){
          message('\nNo box plot of differential peaks from', m, '\n')
        }
      )

      tryCatch(
        {
          profile_colors <- list()
          for (i in 1:length(unique(dbObj.analyzed$samples$Condition))) {
            profile_colors[[unique(dbObj.analyzed$samples$Condition)[i]]] <- c('white', discrete_colours[[i]])
          }
          dba.plotProfile(dbObj.analyzed, merge=c(DBA_REPLICATE), normalize=TRUE)

          dba.plotProfile(profiles, matrices_color=profile_colors, 
                          all_color_scales_equal=FALSE, 
                          decreasing=FALSE, 
                          #group_anno_color=unlist(unname(conditions_colour_code))
                          )
        }, error=function(e){
          message('\nNo profile plot of differential peaks from', m, '\n')
        }
      )

    }
    

  },
  error = function(e) {
    message("\nError occurred: ", conditionMessage(e), "\n")
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)


# ========= Generate DE Analysis Reports =========
tryCatch(
  {
    # Output log to file
    con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
    sink(con, split = FALSE)                 # normal output
    #sink(con, type = "message", split = TRUE)  # messages

    reports <- list()
    reports[["DESeq2"]] <- dba.report(dbObj.analyzed, method=DBA_DESEQ2, contrast=1, th=1)
    reports[["edgeR"]] <- dba.report(dbObj.analyzed, method=DBA_EDGER, contrast=1, th=1)

  },
  error = function(e) {
    message("\nError occurred: ", conditionMessage(e), "\n")
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
)

for (report in names(reports)){
  output_prefix <- change_dirs(result_dir, report, '')
  print("\n\n---------- ", report, " ----------\n\n")
  tryCatch(
    {
      # Output log to file
      con <- file(paste(opt$result_dir, str_replace(opt$result_dir, "/", "_log.txt") , sep=''), open = "at")
      sink(con, split = FALSE)                 # normal output
      #sink(con, type = "message", split = TRUE)  # messages

      pdf(paste(output_prefix, report, '_affinity-analysis_annotated-sites.pdf', sep=""),
          width  = figs$width,
          height = figs$height,
          pointsize = figs$pointsize)  # or 9–10)
      
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
        cat("\nNo DB sites identified by", report, "at a significance threshold of", opt$fdr, "skipping further analysis...\n")
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
      cat('\n\t', dim(gained)[[1]], 'DB peaks in', dbObj.contrast$contrasts[[1]]$name2, '(log2Fold-change > ', opt$lfc, '\n')
      cat('\n\t', dim(lost)[[1]], 'DB peaks in', dbObj.contrast$contrasts[[1]]$name1, '(log2Fold-change < -', opt$lfc, '\n')
      
      
      # Write DE result to bed files
      if (dim(gained)[1] > 0){
        output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name1)
        write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
      }
      if (dim(lost)[1] > 0){
        output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
        write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.bed', sep=''), sep="\t", quote=F, row.names=F, col.names=F)
      }
      

      # Plot profile heatmaps for all significant (FDR < opt$fdr) sites for each method
      output_prefix <- change_dirs(result_dir, report, '')
      tryCatch(
        {
          profile_colors <- list()
          for (i in 1:length(unique(dbObj.analyzed$samples$Condition))) {
            profile_colors[[unique(dbObj.analyzed$samples$Condition)[i]]] <- c('white', discrete_colours[[i]])
          }
          #profile_colors <- rev(profile_colors)
          # Plots all significant sites among both conditions (will include signals)
          dba.plotProfile(dbObj.analyzed, merge=c(DBA_REPLICATE), normalize=TRUE, 
                          sites=GRanges(out %>% dplyr::filter(FDR < opt$fdr))
          )
          dba.plotProfile(profiles_significant, matrices_color=profile_colors, all_color_scales_equal=FALSE, 
                          decreasing=FALSE, 
                          #group_anno_color=unlist(unname(conditions_colour_code))
          )
        }, error=function(e){
          message("\nNo profile plot for ", report, "\n", e)
        }
      )
      
      gained <- out %>% 
        dplyr::filter(FDR < opt$fdr & Fold > opt$lfc)
      lost <- out %>% 
        dplyr::filter(FDR < opt$fdr & Fold < (0 - opt$lfc))
      
      # Write complete DB result to files
      if (dim(gained)[1] > 0){
        output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
        gained <- as.data.frame(annotatePeak(GRanges(gained), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                            level=opt$annotation_level,
                                            tssRegion=c(-3000, 3000))@anno)
        gained <- gained[order(gained$Fold, decreasing=TRUE), ]
        write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name2, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
      }
      if (dim(lost)[1] > 0){
        output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name1)
        lost <- as.data.frame(annotatePeak(GRanges(lost), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                          level=opt$annotation_level,
                                          tssRegion=c(-3000, 3000))@anno)
        lost <- lost[order(lost$Fold, decreasing=FALSE), ]
        write.table(lost, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
      }
      
      # Use absolute fold-change for plotting magnitude
      gained$Fold <- abs(gained$Fold)
      lost$Fold <- abs(lost$Fold)
      
      gpeaks <- GenomicRanges::GRangesList(Gained=gained, Lost=lost)
      names(gpeaks) <- c(dbObj.contrast$contrasts[[1]]$name2, dbObj.contrast$contrasts[[1]]$name1)
      
      # Plot peaks gained/lost over genome between conditions
      output_prefix <- change_dirs(result_dir, report, '')
      #tryCatch(
      #  {
      #    plt <- covplot(gpeaks, title=paste("Peaks over Genome", report, sep=' - '), weightCol='Fold') + 
      #      scale_color_manual(values=rev(c(colours[1:length(unique_peaks)]))) + 
      #      scale_fill_manual(values=rev(c(colours[1:length(unique_peaks)])))
      #    invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_significant_merged_peaks.png', sep=''), plot=plt, dpi=320)))
      #    plt <- plt + facet_grid(chr ~ .id)
      #    invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_significant_peaks.png', sep=''), plot=plt, dpi=320)))
      #  }, error=function(e){
      #    message("No figure\n", e)
      #  }
      #)
      #invisible(capture.output(gc()))
      
      # ========= Get Annotations =========
      peakAnnoList <- list()
      for (p in names(gpeaks)){
        cat("\nObtaining annotations for ", p, "\n")
        output_prefix <- change_dirs(result_dir, report, p)
        
        if (length(gpeaks[[p]]) == 0){
          cat("\nNo", p, "peaks to annotate...\n")
          next
        }

        # Set plot colors
        colour <- conditions_colour_code[[p]]
        heat_colour <- conditions_colour_code[[p]]
        
        if (length(gpeaks[[p]]) > 0){
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
                                        organism=anno_ref$keggOrg,
                                        universe = bg_entrez
              ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
              if (!class(compKEGG) == 'compareClusterResult'){
                cat("\nNo results.\n")
                #next
              }else{
                cat('\n', dim(compKEGG@compareClusterResult)[1], 'results\n')
              }
              if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
                # Map EntrezIDs to gene SYMBOL
                compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
                myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
                for (i in 1:length(myEntrez)){
                  compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                }
                
                plt <- make_anno_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=colour, n=15)
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_KEGG.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
                
                # Write annotations to csv
                write.table(as.data.frame(compKEGG), file=paste(output_prefix, report, '_', p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, ylabel="KEGG Category")
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, ylabel="KEGG Category")
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                remove(compKEGG)
                gc()
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
                                        readable=TRUE,
                                        universe = bg_genes
                ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
                if (!class(compGO) == 'compareClusterResult'){
                  cat("\nNo results.\n")
                  next
                }else{
                  cat('\n', dim(compGO@compareClusterResult)[1], 'results\n')
                }
                if (!is.null(compGO)){
                  #compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
                  #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
                  plt <- make_anno_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=colour, n=15)
                  invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
                  
                  # Write annotations to csv
                  write.table(as.data.frame(compGO), file=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                  
                  plt <- make_pheatmapplot(compGO@compareClusterResult, res, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
                  invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap_by_gene.png', sep=''), plot=plt, dpi=320)))
                  remove(plt)
                  plt <- make_pheatmapplot(compGO@compareClusterResult, res, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
                  invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
                  remove(plt)
                  remove(compGO)
                  gc()
                } else{
                  cat("\nNo annotation results\n")
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
          genes_entrez <- list()
          genes_entrez[[p]] <- anno@anno$geneId
          names(genes_entrez) = sub("_", "\n", names(genes_entrez))
          for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
            cat("\nGetting ", report, ' ', p, ' DAVID ', annotation_type, ' annotations...\n')  
            # DAVID Annotation
            tryCatch(
              {
                cat("\nDAVID - ", annotation_type, "\n")
                entrez <- unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))])
                
                py_run_string("import pandas as pd")
                py_run_string("import sys")
                py_run_string("from suds.client import Client")
                
                # create a service client using the wsdl.
                py_run_string("client = Client('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl')")
                py_run_string("client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')")
                
                #authenticate user email
                py_run_string("client.service.authenticate(r.opt['david_user'])")
                
                # Read input gene list file, convert ids to a comma-delimited string and upload the list to DAVID
                py_run_string("inputIds = ','.join(r.entrez)")
                py_run_string("client.service.addList(inputIds, 'ENTREZ_GENE_ID', r.p, 0)")
                
                # setCategories
                py_run_string("categorySting = str(client.service.setCategories(r.annotation_type))")
                
                #getChartReport
                py_run_string("thd = r.opt['fdr']")
                py_run_string("ct = r.opt['minGSSize']")
                py_run_string("chartReport = client.service.getChartReport(thd,ct)")
                py_run_string("chartRow = len(chartReport)")
                py_run_string("print ('Total chart records:',chartRow)")
                
                if (py$chartRow > 0){
                  if (annotation_type == "KEGG_PATHWAY"){
                    splitter <- ":"
                  }else{
                    splitter <- "~"
                  }
                  # parse chartReport
                  py_run_string("records = pd.DataFrame()")
                  py_run_string("for simpleChartRecord in chartReport:
                                  df = pd.DataFrame(index=[records.shape[0]], data={
                                      'ID' : simpleChartRecord['termName'].split(r.splitter)[0],
                                      'Category' : simpleChartRecord['categoryName'],
                                      'Description' : simpleChartRecord['termName'].split(r.splitter)[1],
                                      'GeneRatio': str(simpleChartRecord['listHits']) + '/' + str(simpleChartRecord['listTotals']), 
                                      'BgRatio': str(simpleChartRecord['popHits']) + '/' + str(simpleChartRecord['popTotals']), 
                                      'pvalue' : simpleChartRecord['ease'],
                                      'p.adjust' : simpleChartRecord['benjamini'],
                                      'FDR' : simpleChartRecord['afdr'],
                                      'geneID' : simpleChartRecord['geneIds'].replace(', ', '/'),
                                      'Count' : simpleChartRecord['listHits'],
                                      'foldEnrichment' : simpleChartRecord['foldEnrichment'],
                                      'id' : simpleChartRecord['id']
                                      })
                                  records = pd.concat([records, df])")
                  #py_run_string("records.reset_index(drop=True, inplace=True)")
                  compD <- py$records
                  
                  qobj <- tryCatch(qvalue(p=as.numeric(compD$pvalue), lambda=opt$pvalue, pi0.method="bootstrap"),
                                  error=function(e) NULL)
                  if (class(qobj) == "qvalue") {
                    qvalues <- qobj$qvalues
                  } else {
                    qvalues <- NA
                  }
                  compD$qvalue <- qvalues
                  
                  compDAVID <- new("enrichResult",
                                  result         = compD,
                                  pvalueCutoff   = opt$fdr,
                                  pAdjustMethod  = "BH",
                                  organism       = opt$assembly,
                                  ontology       = annotation_type,
                                  gene           = entrez,
                                  keytype        = "ENTREZ_GENE_ID")
                  rm(compD)
                  
                  if (!class(compDAVID) == 'enrichResult'){
                    cat("\nNo results.\n")
                    next
                  }else{
                    cat('\n', dim(compDAVID@result)[1], 'results\n')
                  }
                  # old deprecated function
                  #compDAVID <- enrichDAVID(
                  #  unname(genes_entrez[[n]][!is.na(unname(genes_entrez[[n]]))]),
                  #  idType = "ENTREZ_GENE_ID",
                  #  minGSSize = opt$minGSSize,
                  #  maxGSSize = opt$maxGSSize,
                  #  annotation = annotation_type,
                  #  pvalueCutoff = opt$pvalue,
                  #  pAdjustMethod = "BH",
                  #  #species = NA,
                  #  david.user=opt$david_user
                  #)
                  
                  if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
                    # Map EntrezIDs to gene SYMBOL
                    compDAVID@result$SYMBOL <- compDAVID@result$geneID
                    myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
                    for (i in 1:length(myEntrez)){
                      compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                    }
                    
                    plt <- make_anno_dotplot(compDAVID@result, title=paste('DAVID - ', p, sep=""), ylabel=paste(annotation_type,"Category", sep=' '), colour=colour, n=15)
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
                    remove(compDAVID)
                    gc()
                  } else{
                    cat("\nNo annotation results\n")
                    remove(compDAVID)
                    gc()
                  }
                }
              },error = function(e)
              {
                message(e)
                remove(compDAVID)
                gc()
              }
            )
          }
        }else{
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
        gained <- gained[order(gained$Fold.DESeq2, decreasing=TRUE), ]
        write.table(gained, file=paste(output_prefix, 'analyzed_report_', report, '_', dbObj.contrast$contrasts[[1]]$name1, '.tsv', sep=''), sep="\t", quote=F, row.names=F)
      }
      if (dim(lost)[1] > 0){
        output_prefix <- change_dirs(result_dir, report, dbObj.contrast$contrasts[[1]]$name2)
        lost <- as.data.frame(annotatePeak(GRanges(lost), TxDb=anno_ref$txdb, annoDb=anno_ref$annoDb,
                                          level=opt$annotation_level,
                                          tssRegion=c(-3000, 3000))@anno)
        lost <- lost[order(lost$Fold.DESeq2, decreasing=FALSE), ]
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
                                        organism=anno_ref$keggOrg,
                                        universe = bg_entrez
              ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
              if (!class(compKEGG) == 'compareClusterResult'){
                cat("\nNo results.\n")
                #next
              }else{
                cat('\n', dim(compKEGG@compareClusterResult)[1], 'results\n')
              }
              if ((!is.null(compKEGG)) & (dim(compKEGG@compareClusterResult)[1] > 0)){
                # Map EntrezIDs to gene SYMBOL
                compKEGG@compareClusterResult$SYMBOL <- compKEGG@compareClusterResult$geneID
                myEntrez <- lapply(compKEGG@compareClusterResult$geneID, strsplit, '/')
                for (i in 1:length(myEntrez)){
                  compKEGG@compareClusterResult$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                }
                
                plt <- make_anno_dotplot(compKEGG@compareClusterResult, title=paste('KEGG - ', p, sep=""), ylabel="KEGG Category", colour=colour, n=15)
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_KEGG.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
                
                # Write annotations to csv
                write.table(as.data.frame(compKEGG), file=paste(output_prefix, report, '_', p, '_annotated_KEGG.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                
                plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, ylabel="KEGG Category")
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap_bygene.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                plt <- make_pheatmapplot(compKEGG@compareClusterResult, res, anno_type="KEGG", assembly=opt$assembly, title=paste('KEGG - ', p, sep=""), heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, ylabel="KEGG Category")
                invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_KEGG_annotation_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
                remove(plt)
                remove(compKEGG)
                gc()
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
                                        readable=TRUE,
                                        universe = bg_genes
                ) # Check https://www.genome.jp/kegg/catalog/org_list.html for organism hsa=human mmu=mouse
                if (!class(compGO) == 'compareClusterResult'){
                  cat("\nNo results.\n")
                  next
                }else{
                  cat('\n', dim(compGO@compareClusterResult)[1], 'results\n')
                }
                if ((!is.null(compGO)) & (dim(compGO@compareClusterResult)[1] > 0)){
                  #compGO@compareClusterResult$ONTOLOGY <- go2ont(compGO@compareClusterResult$ID)$Ontology
                  #plt <- dotplot(compGO, showCategory = 10, title = "GO Pathway Enrichment Analysis")
                  plt <- make_anno_dotplot(compGO@compareClusterResult, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term", colour=colour, n=15)
                  invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.png', sep=''), plot=plt, dpi=320, width=10, units='in')))
                  
                  # Write annotations to csv
                  write.table(as.data.frame(compGO), file=paste(output_prefix, report, '_', p, '_annotated_GO-', ont, '.tsv', sep=''), sep="\t", quote=F, row.names=F, col.names=T)
                  
                  plt <- make_pheatmapplot(compGO@compareClusterResult, res, assembly=opt$assembly, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=TRUE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
                  invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap_by_gene.png', sep=''), plot=plt, dpi=320)))
                  remove(plt)
                  plt <- make_pheatmapplot(compGO@compareClusterResult, res, assembly=opt$assembly, heat_colour = heat_colour, num_terms=25, num_genes=50, lfc=round(max(res$log2FoldChange)), dendro=TRUE, sort_genes=FALSE, title=paste("GO (", ont, ") - ", p, sep=""), ylabel="GO Term")
                  invisible(capture.output(ggsave(filename=paste(output_prefix, report, '_GO-', ont, '_', p, '_pheatmap.png', sep=''), plot=plt, dpi=320)))
                  remove(plt)
                  remove(compGO)
                  gc()
                } else{
                  cat("\nNo annotation results\n")
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
          
          genes_entrez <- list()
          genes_entrez[[p]] <- anno@anno$geneId
          names(genes_entrez) = sub("_", "\n", names(genes_entrez))
          for (annotation_type in c("GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY")){
            cat("\nGetting ", report, ' ', p, ' DAVID ', annotation_type, ' annotations...\n')  
            # DAVID Annotation
            tryCatch(
              {
                cat("\nDAVID - ", annotation_type, "\n")
                entrez <- unname(genes_entrez[[p]][!is.na(unname(genes_entrez[[p]]))])
                
                py_run_string("import pandas as pd")
                py_run_string("import sys")
                py_run_string("from suds.client import Client")
                
                # create a service client using the wsdl.
                py_run_string("client = Client('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService?wsdl')")
                py_run_string("client.wsdl.services[0].setlocation('https://davidbioinformatics.nih.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/')")
                
                #authenticate user email
                py_run_string("client.service.authenticate(r.opt['david_user'])")
                
                # Read input gene list file, convert ids to a comma-delimited string and upload the list to DAVID
                py_run_string("inputIds = ','.join(r.entrez)")
                py_run_string("client.service.addList(inputIds, 'ENTREZ_GENE_ID', r.p, 0)")
                
                # setCategories
                py_run_string("categorySting = str(client.service.setCategories(r.annotation_type))")
                
                #getChartReport
                py_run_string("thd = r.opt['fdr']")
                py_run_string("ct = r.opt['minGSSize']")
                py_run_string("chartReport = client.service.getChartReport(thd,ct)")
                py_run_string("chartRow = len(chartReport)")
                py_run_string("print ('Total chart records:',chartRow)")
                
                if (py$chartRow > 0){
                  if (annotation_type == "KEGG_PATHWAY"){
                    splitter <- ":"
                  }else{
                    splitter <- "~"
                  }
                  # parse chartReport
                  py_run_string("records = pd.DataFrame()")
                  py_run_string("for simpleChartRecord in chartReport:
                              df = pd.DataFrame(index=[records.shape[0]], data={
                                  'ID' : simpleChartRecord['termName'].split(r.splitter)[0],
                                  'Category' : simpleChartRecord['categoryName'],
                                  'Description' : simpleChartRecord['termName'].split(r.splitter)[1],
                                  'GeneRatio': str(simpleChartRecord['listHits']) + '/' + str(simpleChartRecord['listTotals']), 
                                  'BgRatio': str(simpleChartRecord['popHits']) + '/' + str(simpleChartRecord['popTotals']), 
                                  'pvalue' : simpleChartRecord['ease'],
                                  'p.adjust' : simpleChartRecord['benjamini'],
                                  'FDR' : simpleChartRecord['afdr'],
                                  'geneID' : simpleChartRecord['geneIds'].replace(', ', '/'),
                                  'Count' : simpleChartRecord['listHits'],
                                  'foldEnrichment' : simpleChartRecord['foldEnrichment'],
                                  'id' : simpleChartRecord['id']
                                  })
                              records = pd.concat([records, df])")
                  #py_run_string("records.reset_index(drop=True, inplace=True)")
                  compD <- py$records
                  
                  qobj <- tryCatch(qvalue(p=as.numeric(compD$pvalue), lambda=opt$pvalue, pi0.method="bootstrap"),
                                  error=function(e) NULL)
                  if (class(qobj) == "qvalue") {
                    qvalues <- qobj$qvalues
                  } else {
                    qvalues <- NA
                  }
                  compD$qvalue <- qvalues
                  
                  compDAVID <- new("enrichResult",
                                  result         = compD,
                                  pvalueCutoff   = opt$fdr,
                                  pAdjustMethod  = "BH",
                                  organism       = opt$assembly,
                                  ontology       = annotation_type,
                                  gene           = entrez,
                                  keytype        = "ENTREZ_GENE_ID")
                  rm(compD)
                  
                  if (!class(compDAVID) == 'enrichResult'){
                    cat("\nNo results.\n")
                    next
                  }else{
                    cat('\n', dim(compDAVID@result)[1], 'results\n')
                  }
                  # old deprecated function
                  #compDAVID <- enrichDAVID(
                  #  unname(genes_entrez[[n]][!is.na(unname(genes_entrez[[n]]))]),
                  #  idType = "ENTREZ_GENE_ID",
                  #  minGSSize = opt$minGSSize,
                  #  maxGSSize = opt$maxGSSize,
                  #  annotation = annotation_type,
                  #  pvalueCutoff = opt$pvalue,
                  #  pAdjustMethod = "BH",
                  #  #species = NA,
                  #  david.user=opt$david_user
                  #)
                  
                  if ((!is.null(compDAVID)) & (dim(compDAVID@result)[1] > 0)){
                    # Map EntrezIDs to gene SYMBOL
                    compDAVID@result$SYMBOL <- compDAVID@result$geneID
                    myEntrez <- lapply(compDAVID@result$geneID, strsplit, '/')
                    for (i in 1:length(myEntrez)){
                      compDAVID@result$SYMBOL[i] <- paste(plyr::mapvalues(myEntrez[[i]][[1]], mapper$geneId, mapper$SYMBOL, warn_missing = FALSE), collapse='/')
                    }
                    
                    plt <- make_anno_dotplot(compDAVID@result, title=paste('DAVID - ', p, sep=""), ylabel=paste(annotation_type,"Category", sep=' '), colour=colour, n=15)
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
                    remove(compDAVID)
                    gc()
                  } else{
                    cat("\nNo annotation results\n")
                    remove(compDAVID)
                    gc()
                  }
                }
                
              },error = function(e)
              {
                message(e)
                remove(compDAVID)
                gc()
              }
            )
          }
            
        }else{
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
      
  },
  error = function(e) {
    message("\nError occurred: ", conditionMessage(e), "\n")
  }, 
  finally = {
    if (!is.null(dev.list())) dev.off()
    invisible(capture.output(gc()))
    #sink()                  # stop messages
    sink()                  # stop normal output
    close(con)
  }
  
)