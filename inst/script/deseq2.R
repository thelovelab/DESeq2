# A command-line interface to DESeq2 for use with Galaxy
# written by Bjoern Gruening and modified by Michael Love 2015.01.11
#
# one of these arguments is required:
#
#   'factors' a JSON list object from Galaxy
#
#   'sample_table' is a sample table as described in ?DESeqDataSetFromHTSeq
#   with columns: sample name, filename, then factors (variables)
#
# the output file has columns:
# 
#   baseMean (mean normalized count)
#   log2FoldChange (by default a moderated LFC estimate)
#   lfcSE (the standard error)
#   stat (the Wald statistic)
#   pvalue (p-value from comparison of Wald statistic to a standard Normal)
#   padj (adjusted p-value, Benjamini Hochberg correction on genes which pass the mean count filter)
# 
# the first variable in 'factors' and first column in 'sample_table' will be the primary factor.
# the levels of the primary factor are used in the order of appearance in factors or in sample_table.
#
# by default, levels in the order A,B,C produces a single comparison of B vs A, to a single file 'outfile'
#
# for the 'many_contrasts' flag, levels in the order A,B,C produces comparisons C vs A, B vs A, C vs B,
# to a number of files using the 'outfile' prefix: 'outfile.condition_C_vs_A' etc.
# all plots will still be sent to a single PDF, named by the arg 'plots', with extra pages.
#
# fit_type is an integer valued argument, with the options from ?estimateDisperions
#   1 "parametric"
#   2 "local"
#   3 "mean"

# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("getopt")
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)
args <- commandArgs(trailingOnly = TRUE)

# get options, using the spec as defined by the enclosed list.
# we read the options from the default: commandArgs(TRUE).
spec <- matrix(c(
  "quiet", "q", 0, "logical",
  "help", "h", 0, "logical",
  "outfile", "o", 1, "character",
  "factors", "f", 1, "character",
  "plots" , "p", 1, "character",
  "sample_table", "s", 1, "character",
  "fit_type", "t", 1, "integer",
  "many_contrasts", "m", 0, "logical",
  "outlier_replace_off" , "a", 0, "logical",
  "outlier_filter_off" , "b", 0, "logical",
  "auto_mean_filter_off", "c", 0, "logical",
  "beta_prior_off", "d", 0, "logical"),
  byrow=TRUE, ncol=4)
opt <- getopt(spec)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

# enforce the following required arguments
if (is.null(opt$outfile)) {
  cat("'outfile' is required\n")
  q(status=1)
}
if (is.null(opt$sample_table) & is.null(opt$factors)) {
  cat("'factors' or 'sample_table' is required\n")
  q(status=1)
}

verbose <- if (is.null(opt$quiet)) {
  TRUE
} else {
  FALSE
}

suppressPackageStartupMessages({
  library("DESeq2")
  library("RColorBrewer")
  library("gplots")
})

# build or read sample table

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

# switch on if 'factors' was provided:
if (!is.null(opt$factors)) {
  library("rjson")
  parser <- newJSONParser()
  parser$addData(opt$factors)
  factorList <- parser$getObject()
  factors <- sapply(factorList, function(x) x[[1]])
  primaryFactor <- factors[1]
  filenamesIn <- unname(unlist(factorList[[1]][[2]]))
  sampleTable <- data.frame(sample=basename(filenamesIn), filename=filenamesIn, row.names=filenamesIn)
  for (factor in factorList) {
    factorName <- trim(factor[[1]])
    sampleTable[[factorName]] <- character(nrow(sampleTable))
    lvls <- sapply(factor[[2]], function(x) names(x))
    for (i in seq_along(factor[[2]])) {
      files <- factor[[2]][[i]][[1]]
      sampleTable[files,factorName] <- trim(lvls[i])
    }
    sampleTable[[factorName]] <- factor(sampleTable[[factorName]], levels=lvls)
  }
  rownames(sampleTable) <- sampleTable$sample
} else {
  # read the sample_table argument
  # this table is described in ?DESeqDataSet
  # one column for the sample name, one for the filename, and
  # the remaining columns for factors in the analysis
  sampleTable <- read.delim(opt$sample_table)
  factors <- colnames(sampleTable)[-c(1:2)]
  for (factor in factors) {
    lvls <- unique(as.character(sampleTable[[factor]]))
    sampleTable[[factor]] <- factor(sampleTable[[factor]], levels=lvls)
  }
}

primaryFactor <- factors[1]
designFormula <- as.formula(paste("~", paste(rev(factors), collapse=" + ")))

if (verbose) {
  cat("DESeq2 run information\n\n")
  cat("sample table:\n")
  print(sampleTable[,-c(1:2),drop=FALSE])
  cat("\ndesign formula:\n")
  print(designFormula)
  cat("\n\n")
}

# these are plots which are made once for each analysis
generateGenericPlots <- function(dds, factors) {
  rld <- rlog(dds)
  print(plotPCA(rld, intgroup=rev(factors)))
  # need meaningful labels, because from Galaxy, sample names are random
  labs <- paste0(seq_len(ncol(dds)), ": ", do.call(paste, as.list(colData(dds)[factors])))
  dat <- assay(rld)
  colnames(dat) <- labs
  distsRL <- dist(t(dat))
  mat <- as.matrix(distsRL)
  hc <- hclust(distsRL)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol),
            main="Sample-to-sample distances", margin=c(13,13))
  plotDispEsts(dds, main="Dispersion estimates")
}

# these are plots which can be made for each comparison, e.g.
# once for C vs A and once for B vs A
generateSpecificPlots <- function(res, threshold, title_suffix) {
  use <- res$baseMean > threshold
  if (sum(!use) == 0) {
    h <- hist(res$pvalue, breaks=0:50/50, plot=FALSE)
    barplot(height = h$counts,
            col = "powderblue", space = 0, xlab="p-values", ylab="frequency",
            main=paste("Histogram of p-values for",title_suffix))
    text(x = c(0, length(h1$counts)), y = 0, label=paste(c(0,1)), adj=c(0.5,1.7), xpd=NA)
  } else {
    h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
    h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
    colori <- c("filtered (low count)"="khaki", "not filtered"="powderblue")
    barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
            col = colori, space = 0, xlab="p-values", ylab="frequency",
            main=paste("Histogram of p-values for",title_suffix))
    text(x = c(0, length(h1$counts)), y = 0, label=paste(c(0,1)), adj=c(0.5,1.7), xpd=NA)
    legend("topright", fill=rev(colori), legend=rev(names(colori)), bg="white")
  }
    plotMA(res, main= paste("MA-plot for",title_suffix), ylim=range(res$log2FoldChange, na.rm=TRUE))
}

if (verbose) {
  cat(paste("primary factor:",primaryFactor,"\n"))
  if (length(factors) > 1) {
    cat(paste("other factors in design:",paste(factors[-length(factors)],collapse=","),"\n"))
  }
  cat("\n---------------------\n")
}

# if JSON input from Galaxy, path is absolute
# otherwise, from sample_table, assume it is relative
dir <- if (is.null(opt$factors)) {
  "."
} else {
  ""
}

# construct the object
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = dir,
                                  design =  designFormula)

if (verbose) cat(paste(ncol(dds), "samples with counts over", nrow(dds), "genes\n"))

# optional outlier behavior
if (is.null(opt$outlier_replace_off)) {
  minRep <- 7
} else {
  minRep <- Inf
  if (verbose) cat("outlier replacement off\n")
}
if (is.null(opt$outlier_filter_off)) {
  cooksCutoff <- TRUE
} else {  
  cooksCutoff <- FALSE
  if (verbose) cat("outlier filtering off\n")
}

# optional automatic mean filtering
if (is.null(opt$auto_mean_filter_off)) {
  independentFiltering <- TRUE
} else {
  independentFiltering <- FALSE
  if (verbose) cat("automatic filtering on the mean off\n")
}

# shrinkage of LFCs
if (is.null(opt$beta_prior_off)) {
  betaPrior <- TRUE
} else {
  betaPrior <- FALSE
  if (verbose) cat("beta prior off\n")
}

# dispersion fit type
if (is.null(opt$fit_type)) {
  fitType <- "parametric"
} else {
  fitType <- c("parametric","local","mean")[opt$fit_type]
}

if (verbose) cat(paste("using disperion fit type:",fitType,"\n"))

# run the analysis
dds <- DESeq(dds, fitType=fitType, betaPrior=betaPrior, minReplicatesForReplace=minRep)

# create the generic plots and leave the device open
if (!is.null(opt$plots)) {
  if (verbose) cat("creating plots\n")
  pdf(opt$plots)
  generateGenericPlots(dds, factors)
}

n <- nlevels(colData(dds)[[primaryFactor]])
allLevels <- levels(colData(dds)[[primaryFactor]])

if (is.null(opt$many_contrasts)) {
  # only contrast the first and second level of the primary factor
  ref <- allLevels[1]
  lvl <- allLevels[2]
  res <- results(dds, contrast=c(primaryFactor, lvl, ref),
                 cooksCutoff=cooksCutoff,
                 independentFiltering=independentFiltering)
  if (verbose) {
    cat("summary of results\n")
    cat(paste0(primaryFactor,": ",lvl," vs ",ref,"\n"))
    print(summary(res))
  }
  resSorted <- res[order(res$padj),]
  outDF <- as.data.frame(resSorted)
  outDF$geneID <- rownames(outDF)
  outDF <- outDF[,c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  filename <- opt$outfile
  write.table(outDF, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  if (independentFiltering) {
    threshold <- unname(attr(res, "filterThreshold"))
  } else {
    threshold <- 0
  }
  title_suffix <- paste0(primaryFactor,": ",lvl," vs ",ref)
  if (!is.null(opt$plots)) {
    generateSpecificPlots(res, threshold, title_suffix)
  }
} else {
  # rotate through the possible contrasts of the primary factor
  # write out a sorted table of results with the contrast as a suffix
  # add contrast specific plots to the device
  for (i in seq_len(n-1)) {
    ref <- allLevels[i]
    contrastLevels <- allLevels[(i+1):n]
    for (lvl in contrastLevels) {
      res <- results(dds, contrast=c(primaryFactor, lvl, ref),
                     cooksCutoff=cooksCutoff,
                     independentFiltering=independentFiltering)
      resSorted <- res[order(res$padj),]
      outDF <- as.data.frame(resSorted)
      outDF$geneID <- rownames(outDF)
      outDF <- outDF[,c("geneID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
      filename <- paste0(opt$outfile,".",primaryFactor,"_",lvl,"_vs_",ref)
      write.table(outDF, file=filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
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
  }
}

# close the plot device
if (!is.null(opt$plots)) {
  cat("closing plot device\n")
  dev.off()
}

cat("Session information:\n\n")

sessionInfo()

