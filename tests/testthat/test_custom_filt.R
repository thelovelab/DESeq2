# try a custom filter function
set.seed(1)
dds <- makeExampleDESeqDataSet(n=200, m=4, betaSD=rep(c(0,2),c(150,50)))
dds <- DESeq(dds)
res <- results(dds)

filter <- mcols(dds)$baseMean
test <- res$pvalue
theta <- seq(mean(filter == 0), 1, length=20)
method <- "BH"
alpha <- 0.1

customFilt <- function(filter, test, theta, method, alpha) {
  cutoff <- quantile(filter, theta)
  numRej <- sapply(cutoff, function(x) sum(p.adjust(test[filter > x]) < alpha, na.rm=TRUE))
  threshold <- theta[which(numRej == max(numRej))[1]]
  padj <- numeric(length(test))
  padj <- NA
  idx <- filter > quantile(filter, threshold)
  padj[idx] <- p.adjust(test[idx], method=method)
  return(padj)
}

resCustom <- results(dds, filterFun=customFilt)
#plot(res$padj, resCustom$padj);abline(0,1)
