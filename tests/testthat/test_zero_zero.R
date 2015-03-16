# test comparison of two groups with all zeros
dds <- makeExampleDESeqDataSet(m=8, n=100, sizeFactors=c(1,1,.5,.5,1,1,2,2))
dds$condition <- factor(rep(c("A","B","C","D"),each=2))
counts(dds)[1,] <- c(100L,110L,0L,0L,100L,110L,0L,0L)
counts(dds)[2,] <- rep(0L, 8)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","D","B"))[1,]
expect_equal(res$log2FoldChange, 0)
res <- results(dds, contrast=c(0,0,-1,0,1))[1,]
expect_equal(res$log2FoldChange, 0)
res <- results(dds,c(0,0,0,0,1))[1,]
expect_true(res$log2FoldChange != 0)
# if all samples have 0, should be NA
res <- results(dds, contrast=c("condition","D","B"))[2,]
expect_true(is.na(res$log2FoldChange))
res <- results(dds, contrast=c(0,0,-1,0,1))[2,]
expect_true(is.na(res$log2FoldChange))

