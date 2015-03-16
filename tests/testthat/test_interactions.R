dds <- makeExampleDESeqDataSet(n=100,m=8)
colData(dds)$group <- factor(rep(c("x","y"),times=ncol(dds)/2))
design(dds) <- ~ condition + condition:group
dds <- DESeq(dds)
design(dds) <- ~ condition + group + condition:group
dds <- DESeq(dds)

# second order interactions error
dds <- makeExampleDESeqDataSet(n=100,m=16)
dds$x <- factor(rep(1:2,each=8))
dds$y <- factor(rep(rep(1:2,2),each=4))
dds$z <- factor(rep(rep(1:2,4),each=2))
design(dds) <- ~ x*y*z
expect_error(DESeq(dds))

