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
beta0 <- seq(from=1,to=16,length=100)
idx <- rep(rep(c(TRUE,FALSE),c(1,9)),10)
set.seed(1)
#par(mfrow=c(2,3))
for (disp0 in c(.01,.5)) {
  for (m in c(10,20,80)) {
    dds <- makeExampleDESeqDataSet(n=100, m=m, interceptMean=beta0, interceptSD=0,
                                   dispMeanRel=function(x) 4/x + disp0)
    counts(dds)[idx,1] <- as.integer(1000 * 2^beta0[idx])
    dds <- DESeq(dds, minReplicatesForReplace=Inf, quiet=TRUE)
    res <- results(dds)
    cutoff <- qf(.99, 2, m-2)
    outlierCooks <- assays(dds)[["cooks"]][idx,1] > cutoff
    maxOtherCooks <- apply(assays(dds)[["cooks"]][idx,-1], 1, max) < cutoff
    expect_true(all(is.na(res$pvalue[idx])))
    expect_true(all(outlierCooks))
    expect_true(all(maxOtherCooks))
    #col <- rep("black", 100)
    #col[idx] <- ifelse(outlierCooks, ifelse(maxOtherCooks, "blue", "red"), "purple")
    #plot(assays(dds)[["cooks"]][,1], col=col, log="y",
    #     main=paste(m,"-",disp0), ylab="cooks");abline(h=qf(.99,2,m-2))
  }
}
