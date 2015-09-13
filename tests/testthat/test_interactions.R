dds <- makeExampleDESeqDataSet(n=100,m=8)
colData(dds)$group <- factor(rep(c("X","Y"),times=ncol(dds)/2))
design(dds) <- ~ condition + group + condition:group
dds <- DESeq(dds)
expect_equal(resultsNames(dds)[4], "conditionB.groupY")
# interactions error
expect_error(DESeq(dds, betaPrior=TRUE))

