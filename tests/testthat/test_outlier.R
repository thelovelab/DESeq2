# test filtering and replacement
set.seed(1)
dds <- makeExampleDESeqDataSet(n=100, m=12, dispMeanRel = function(x) 4/x + .5)
counts(dds)[1,] <- rep(0L, 12)
counts(dds)[2,] <- c(100000L, rep(10L, 11))
counts(dds)[3,] <- c(100000L, rep(0L, 11))
dds0 <- DESeq(dds, minReplicatesForReplace=Inf)
dds1 <- DESeq(dds, minReplicatesForReplace=6)
pval0 <- results(dds0)[1:3,"pvalue"]
pval <- results(dds1)[1:3,"pvalue"]
LFC0 <- results(dds0)[1:3,"log2FoldChange"]
LFC <- results(dds1)[1:3,"log2FoldChange"]

# filtered
expect_true(all(is.na(pval0)))
# not filtered
expect_true(all(!is.na(pval[2:3])))
# counts still the same
expect_true(all(counts(dds1)==counts(dds)))
# first is NA
expect_true(is.na(LFC[1]))
# replaced, reduced LFC
expect_true(abs(LFC[2]) < abs(LFC0[2]))
# replaced, LFC now zero
expect_true(LFC[3] == 0)
idx <- which(!mcols(dds1)$replace)
# the pvalue for those not replaced is equal
expect_equal(results(dds1)$pvalue[idx], results(dds0)$pvalue[idx])


# check that outlier filtering catches throughout range of mu
#par(mfrow=c(2,2))
beta0 <- seq(from=1,to=16,length=100)
idx <- rep(rep(c(TRUE,FALSE),c(1,9)),10)
set.seed(1)
for (m in c(10,20,60,80)) {
  dds <- makeExampleDESeqDataSet(n=100, m=m, interceptMean=beta0,
                                 dispMeanRel=function(x) 4/x + .5)
  counts(dds)[idx,1] <- as.integer(1000 * 2^beta0[idx])
  dds <- DESeq(dds, minReplicatesForReplace=Inf)
  #plot(assays(dds)[["cooks"]][,1], col=idx + 1, log="y", main=m, ylab="cooks");abline(h=qf(.99,2,m-2))
  res <- results(dds)
  expect_true(all(is.na(res$pvalue[idx])))
}
