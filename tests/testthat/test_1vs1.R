dds <- makeExampleDESeqDataSet(n=100, m=2)
expect_warning({ dds <- DESeq(dds)})
res <- results(dds)
