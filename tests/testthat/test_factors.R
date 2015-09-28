dds <- makeExampleDESeqDataSet(n=100, m=6)
levels(dds$condition) <- c("test-","test+")
expect_error(DESeq(dds))

dds <- makeExampleDESeqDataSet(n=100, m=6)
dds$condition <- factor(rep(letters[1:3], each=2), ordered=TRUE)
expect_error(DESeq(dds))
mm <- model.matrix(~ condition, data=colData(dds))
dds <- DESeq(dds, full=mm) # betaPrior=FALSE
