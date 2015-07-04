dds <- makeExampleDESeqDataSet(n=100, m=6)
levels(dds$condition) <- c("test-","test+")
expect_error(DESeq(dds))
