# one row
set.seed(1)
dds <- makeExampleDESeqDataSet(n=1)
sizeFactors(dds) <- rep(1,ncol(dds))
dispersions(dds) <- .5
dds <- nbinomWaldTest(dds)
res <- results(dds)
dds <- nbinomLRT(dds, reduced=~1)
res <- results(dds)


# only intercept
set.seed(1)
dds <- makeExampleDESeqDataSet(n=100)
design(dds) <- ~ 1
dds <- DESeq(dds)
res <- results(dds)


# metadata insertion
dds <- makeExampleDESeqDataSet(n=50,m=4)

dds2 <- DESeqDataSetFromMatrix( counts(dds), colData(dds), design(dds) )
mcols(dds2)$foo <- paste( "bar", 1:nrow(dds2) )
dds2 <- DESeq(dds2)
results(dds2)
expect_true(class(mcols(mcols(dds2))$type) == "character")

dds3 <- DESeqDataSetFromMatrix( counts(dds), DataFrame(row.names=1:ncol(dds)), ~ 1 )
dds3$test <- 1:ncol(dds3)
dds3 <- estimateSizeFactors(dds3)
expect_true(class(mcols(colData(dds3))$type) == "character")


# underscores
dds <- makeExampleDESeqDataSet(n=50,m=4)
levels(dds$condition) <- c("A_1","B_2")
dds$exp_cond <- dds$condition
design(dds) <- ~ exp_cond
dds <- DESeq(dds)
results(dds)

