# test fast
dds <- makeExampleDESeqDataSet(n=1000, m=20)
rld <- rlog(dds, fast=TRUE)

# expect warning on sparsity and large counts
dds <- makeExampleDESeqDataSet(n=100, m=20)
idx <- sample(ncol(dds), nrow(dds)/2, TRUE)
counts(dds)[cbind(1:(nrow(dds)/2), idx)] <- 10000L
mcols(dds)$dispFit <- .5
expect_warning({ rld <- rlog(dds, blind=FALSE) })
