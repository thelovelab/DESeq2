## A command-line interface to DESeq2 for use with Galaxy
## written by Bjoern Gruening and modified by Michael Love 2015.01.07

# most important:
#
# 'sample_table' is a sample table as described in ?DESeqDataSetFromHTSeq
# with columns: sample name, filename, then factors (variables)
#
# 'factors' is a comma separated list of factors to use in the design
# where the *last* factor is the primary factor for comparisons
#
# 'reference_level' is the level of the primary factor which will be
# compared against, e.g. the control or untreated samples

# Setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "verbose", "v", 0, "logical",
  "help", "h", 0, "logical",
  "outfile_prefix", "o", 1, "character",
  "sample_table", "s", 1, "character",
  "factors", "f", 1, "character",
  "reference_level", "r", 1, "character",
  "plots" , "p", 1, "character",
  "outlier_replace_off" , "a", 0, "logical",
  "outlier_filter_off" , "b", 0, "logical",
  "auto_mean_filter_off", "c", 0, "logical"
  ), byrow=TRUE, ncol=4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# enforce the following required arguments
if (is.null(opt$outfile_prefix)) {
  print("'outfile_prefix' is required")
  q(status=1)
}
if (is.null(opt$sample_table)) {
  print("'sample_table' is required")
  q(status=1)
}
if (is.null(opt$factors)) {
  print("'factors' is required")
  q(status=1)
}

verbose <- if (is.null(opt$verbose)) {
  FALSE
} else {
  TRUE
}

suppressPackageStartupMessages({
  library("DESeq2")
  library("RColorBrewer")
  library("gplots")
})

# these are plots which are made once for each analysis
generateGenericPlots <- function(dds, factors) {
  rld <- rlog(dds)
  print(plotPCA(rld, intgroup=rev(factors)))
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  hc <- hclust(distsRL)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, Rowv=as.dendrogram(hc), trace="none", col = rev(hmcol),
            main="Sample-to-sample distances", margin=c(13,13))
  plotDispEsts(dds, main="Dispersion estimates")
}

# these are plots which are made for each comparison, e.g.
# once for C vs A and once for B vs A
generateSpecificPlots <- function(res, threshold, title_suffix) {
  use <- res$baseMean > threshold
  h1 <- hist(res$pvalue[!use], breaks=50, plot=FALSE)
  h2 <- hist(res$pvalue[use], breaks=50, plot=FALSE)
  colori <- c("filtered (low count)"="khaki", "not filtered"="powderblue")
  barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
          col = colori, space = 0, xlab="p-values", ylab="frequency",
          main=paste("Histogram of p-values for",title_suffix))
  text(x = c(0, length(h1$counts)), y = 0, label=paste(c(0,1)), adj=c(0.5,1.7), xpd=NA)
  legend("topright", fill=rev(colori), legend=rev(names(colori)))
  plotMA(res, main= paste("MA-plot for",title_suffix), ylim=range(res$log2FoldChange, na.rm=TRUE))
}

# trim whitespace on the factors input
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
factors <- unname(sapply(strsplit(opt$factors, ",")[[1]], trim))

# read the sample_table argument
# this table is described in ?DESeqDataSet
# one column for the sample name, one for the filename, and
# the remaining columns for factors in the analysis
sampleTable <- read.delim(opt$sample_table)
designFormula <- as.formula(paste("~", paste(factors, collapse=" + ")))

# make sure the factors listed are columns of sample_table
# and ensure they are interpreted as factors
for (factor in factors) {
  if (!factor %in% colnames(sampleTable)) {
    message(paste(factor,"is not a column of 'sample_table'"))
    q(status=1)
  }
  sampleTable[[factor]] <- factor(sampleTable[[factor]])
}

# the primary factor is the *last* one in the list (to accord with DESeq2, edgeR, limma, etc.)
# the reference level is the *first level alphabetically*, or the one set by argument
primaryFactor <- factors[length(factors)]
if (!is.null(opt$reference_level)) {
  ref <- opt$reference_level
  if (! ref %in% levels(sampleTable[[primaryFactor]])) {
    message(paste(ref,"is not a level of",primaryFactor))
    q(status=1)
  }
  sampleTable[[primaryFactor]] <- relevel(sampleTable[[primaryFactor]], ref=ref)
} else {
  ref <- levels(sampleTable[[primaryFactor]])[1]
}

if (verbose) message(paste("primary factor:",primaryFactor))
if (verbose) message(paste("comparisons against reference level:",ref))
if (length(factors) > 1) {
  if (verbose) message(paste("other factors in design:",paste(factors[-length(factors)],collapse=",")))
}

# construct the object
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = ".",
                                  design =  designFormula)

if (verbose) message(paste(ncol(dds), "samples with counts over", nrow(dds), "genes"))

# optional outlier behavior
if (is.null(opt$outlier_replace_off)) {
  minRep <- 7
} else {
  minRep <- Inf
  if (verbose) message("outlier replacement off")
}
if (is.null(opt$outlier_filter_off)) {
  cooksCutoff <- TRUE
} else {  
  cooksCutoff <- FALSE
  if (verbose) message("outlier filtering off")
}

# optional automatic mean filtering
if (is.null(opt$auto_mean_filter_off)) {
  independentFiltering <- TRUE
} else {
  independentFiltering <- FALSE
  if (verbose) message("automatic filtering on the mean off")
}

# run the analysis
dds <- DESeq(dds, minReplicatesForReplace=minRep, quiet=!verbose)

# create the generic plots and leave the device open
if (!is.null(opt$plots)) {
  if (verbose) message("creating plots")
  pdf(opt$plots)
  generateGenericPlots(dds, factors)
}

n <- nlevels(colData(dds)[[primaryFactor]])
contrastLevels <- levels(colData(dds)[[primaryFactor]])[-1]

# rotate through the n-1 contrasts of the primary factor
# write out a sorted table of results with the contrast as a suffix
# add contrast specific plots to the device
for (lvl in contrastLevels) {
  res <- results(dds, contrast=c(primaryFactor, lvl, ref),
                 cooksCutoff=cooksCutoff,
                 independentFiltering=independentFiltering)
  resSorted <- res[order(res$padj),]
  outDF <- as.data.frame(resSorted)
  outDF$geneID <- rownames(outDF)
  outDF <- outDF[,c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  filename <- paste0(opt$outfile_prefix,".",primaryFactor,"_",lvl,"_vs_",ref)
  write.table(outDF, file=filename, sep="\t", quote=FALSE)
  if (independentFiltering) {
    threshold <- unname(attr(res, "filterThreshold"))
  } else {
    threshold <- 0
  }
  title_suffix <- paste0(primaryFactor,": ",lvl," vs ",ref)
  if (!is.null(opt$plots)) {
    generateSpecificPlots(res, threshold, title_suffix)
  }
}

# close the plot device
if (!is.null(opt$plots)) {
  dev.off()
}

sessionInfo()

