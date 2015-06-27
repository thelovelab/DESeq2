dds <- makeExampleDESeqDataSet(n=100, m=4)
dds$group <- factor(c(1,2,1,2))
design(dds) <- ~ condition
expect_error(DESeq(dds, test="LRT", reduced=~group))
