dds <- makeExampleDESeqDataSet(n=100, m=8)
expect_error(DESeq(dds, test="LRT"))
expect_error(DESeq(dds, test="Wald", full=~condition, reduced=~1))
expect_error(DESeq(dds, full=~1))

m <- model.matrix(~ condition, colData(dds))
expect_error(DESeq(dds, test="LRT", full=m, reduced=~1))
expect_error(DESeq(dds, test="LRT", full=m, reduced=m))
expect_error(DESeq(dds, full=m, betaPrior=TRUE))

design(dds) <- ~ 0 + condition
expect_error(DESeq(dds, betaPrior=TRUE))
