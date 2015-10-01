# expect warning on sparsity and large counts
dds <- makeExampleDESeqDataSet(n=100, m=20)
idx <- sample(ncol(dds), nrow(dds)/2, TRUE)
counts(dds)[cbind(1:(nrow(dds)/2), idx)] <- 10000L
mcols(dds)$dispFit <- .5
expect_warning({ rld <- rlog(dds, blind=FALSE) })

# test rlog basics/errors
dds <- makeExampleDESeqDataSet(n=20, m=4)
colnames(dds) <- NULL
rlog(dds)
head(rlog(assay(dds)))
expect_error(rlog(dds, intercept=rep(1,10)))

mcols(dds)$dispFit <- rep(.5, 20)
rlog(dds, blind=FALSE, intercept=rep(1,20))

expect_error(rlogData(dds))
expect_error(rlogData(dds, intercept=rep(1,10)))

# test normTranform
dds <- makeExampleDESeqDataSet(n=50, m=10)
nt <- normTransform(dds)
plotPCA(nt)
