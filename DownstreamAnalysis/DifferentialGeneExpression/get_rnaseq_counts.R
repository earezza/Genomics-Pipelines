# Script wrapper to generate a counts matrix for RNA-Seq data from bam files
# See https://rdrr.io/bioc/Rsubread/man/featureCounts.html for custom options

# ======= Load Packages =======
suppressWarnings(suppressPackageStartupMessages({
  library(Rsubread)
  library(optparse)
}))

# ======= Get command-line optional arguments =======
option_list = list(
  make_option(c("-b", "--bams"), type="character", default='bams/', help="Path to folder containing all .bam and .bai files for count matrix", metavar="character"),
  make_option(c("-o", "--organism"), type="character", default="mouse", help="Organism to annotate genes (mouse (mm10) or human (hg38) or rat (rn6))", metavar="character"),
  make_option(c("-a", "--assembly"), type="character", default="mm10", help="Assembly to annotate genes (in-built option include mm9, mm10, hg38, hg19), use .gtf file if not in-built", metavar="character"),
  make_option(c("-g", "--gtf"), action="store_true", type="logical", default=FALSE, help="Flag to use a provided .gtf assembly annotation file instead of in-built assembly options", metavar="character"),
  make_option(c("-e", "--ext_annotation"), type="character", default=NULL, help="Path to assembly annotation file if used", metavar="character"),
  make_option(c("-s", "--strand_specific"), type="integer", default=0, help="Obtain counts for strand-specific mappings (0 = not strand specific, 1 = forward-strand only, 2 = reverse-strand only)", metavar="character"),
  make_option(c("-m", "--multi_mappings"), action="store_true", type="logical", default=FALSE, help="Flag to count multi-mapping reads", metavar="character"),
  make_option(c("-d", "--duplicates"), action="store_true", type="logical", default=FALSE, help="Flag to include duplicates", metavar="character"),
  make_option(c("-p", "--paired_end"), action="store_true", type="logical", default=FALSE, help="Flag to indicate if data is paired-end sequences", metavar="character"),
  make_option(c("-r", "--result_dir"), type="character", default="Counts/", help="Directory name for saving output results", metavar="character"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="Number of CPU threads to use when computing counts to speed up run", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cat("Command-line options:\n")
for (i in which(names(opt) != "help")) {
  cat(names(opt)[i], '=', paste(opt)[i], "\n")
}

if (tolower(opt$organism) == "mouse"){
  library(org.Mm.eg.db)
  mapDB <- org.Mm.eg.db
} else if (tolower(opt$organism) == "human"){
  library(org.Hs.eg.db)
  mapDB <- org.Hs.eg.db
} else if (tolower(opt$organism) == "rat"){
  library(org.Rn.eg.db)
  mapDB <- org.Rn.eg.db
} else{
  stop("Invalid choice of organism")
}

if (substr(opt$result_dir, nchar(opt$result_dir), nchar(opt$result_dir)) != "/"){
  opt$result_dir = paste(opt$result_dir, "/", sep="")
}
if (!file.exists(opt$result_dir)) {
  dir.create(opt$result_dir)
}

files <- list.files(path=opt$bams, full.names=TRUE)
bam_files <- files[sapply(files, endsWith, ".bam") == TRUE]

bamcounts <- featureCounts(bam_files, 
                           annot.inbuilt = opt$assembly, 
                           annot.ext = opt$ext_annotation,
                           isGTFAnnotationFile = opt$gtf,
                           countMultiMappingReads = opt$multi_mappings, 
                           ignoreDup = opt$duplicates, 
                           isPairedEnd = opt$paired_end, 
                           nthreads = opt$threads, 
                           verbose = TRUE)

rownames(bamcounts$counts) <- mapIds(mapDB, keys = rownames(bamcounts$counts), column = "SYMBOL", keytype = "ENTREZID")

bamcounts$counts <- bamcounts$counts[!(is.na(rownames(bamcounts$counts))), ]

for (n in names(bamcounts)){
  write.table(bamcounts[[n]], file=paste(opt$result_dir, n, ".csv", sep=""), sep=",", quote=F, col.names=NA)
}
