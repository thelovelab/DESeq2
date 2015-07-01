dds0 <- makeExampleDESeqDataSet(n=100, m=6)
expect_error(levels(dds0$condition) <- c("test-","test+"))

