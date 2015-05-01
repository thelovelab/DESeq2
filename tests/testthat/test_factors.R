dds0 <- makeExampleDESeqDataSet(n=100, m=6)
levels(dds0$condition) <- c("test-","test+")
expect_error({ DESeqDataSetFromMatrix(counts(dds0), colData(dds0), ~ condition) })
