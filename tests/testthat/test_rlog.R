# expect warning on sparsity and large counts
dds <- makeExampleDESeqDataSet(n=100, m=20)
idx <- sample(ncol(dds), nrow(dds)/2, TRUE)
counts(dds)[cbind(1:(nrow(dds)/2), idx)] <- 10000L
mcols(dds)$dispFit <- .5
expect_warning({ rld <- rlog(dds, blind=FALSE) })

# test normTranform
dds <- makeExampleDESeqDataSet(n=50, m=10)
nt <- normTransform(dds)
plotPCA(nt)
