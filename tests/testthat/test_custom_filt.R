context("custom_filt")
test_that("custom filters can be provided to results()", {
  # try a custom filter function
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200, m=4, betaSD=rep(c(0,2),c(150,50)))
  dds <- DESeq(dds)
  res <- results(dds)
  method <- "BH"
  alpha <- 0.1

  customFilt <- function(res, filter, alpha, method) {
    if (missing(filter)) {
      filter <- res$baseMean
    }
    theta <- 0:10/10
    cutoff <- quantile(filter, theta)
    numRej <- sapply(cutoff, function(x) sum(p.adjust(res$pvalue[filter > x]) < alpha, na.rm=TRUE))
    threshold <- theta[which(numRej == max(numRej))[1]]
    res$padj <- numeric(nrow(res))
    idx <- filter > quantile(filter, threshold)
    res$padj[!idx] <- NA
    res$padj[idx] <- p.adjust(res$pvalue[idx], method=method)
    res
  }

  resCustom <- results(dds, filterFun=customFilt)
  plot(res$padj, resCustom$padj);abline(0,1)
})
