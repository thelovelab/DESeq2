coldata <- DataFrame(x=factor(c("A","A","B","B")),
                     name=letters[1:4],
                     ident=factor(rep("A",4)),
                     num=1:4,
                     missinglevels=factor(c("A","A","B","B"), levels=c("A","B","C")),
                     notref=factor(c("control","control","abc","abc")),
                     row.names=1:4)
counts <- matrix(1:16, ncol=4)

expect_message(DESeqDataSet(SummarizedExperiment(list(foo=counts), colData=coldata), ~ x))
expect_error(DESeqDataSetFromMatrix(matrix(c(1:11,-1),ncol=4), coldata, ~ x))
expect_error(DESeqDataSetFromMatrix(matrix(c(1:11,0.5),ncol=4), coldata, ~ x))
expect_error(DESeqDataSetFromMatrix(matrix(rep(0,16),ncol=4), coldata, ~ x))
expect_warning(DESeqDataSetFromMatrix(matrix(rep(1:4,4),ncol=4), coldata, ~ x))
expect_warning(DESeqDataSetFromMatrix(matrix(1:16, ncol=4, dimnames=list(c(1,2,3,3),1:4)), coldata, ~ x))
expect_error(DESeqDataSetFromMatrix(counts, coldata, ~ y))
expect_warning(DESeqDataSetFromMatrix(counts, coldata, ~ name))
expect_error(DESeqDataSetFromMatrix(counts, coldata, ~ ident))
expect_message(DESeqDataSetFromMatrix(counts, coldata, ~ num))
expect_message(DESeqDataSetFromMatrix(counts, coldata, ~ missinglevels))
expect_message(DESeqDataSetFromMatrix(counts, coldata, ~ notref))
expect_error(DESeqDataSetFromMatrix(counts, coldata, ~ident + x), "design contains")

# same colnames but in different order:
expect_error(DESeqDataSetFromMatrix(matrix(1:16, ncol=4, dimnames=list(1:4, 4:1)), coldata, ~ x))

# testing incoming metadata columns
coldata <- DataFrame(x=factor(c("A","A","B","B")))
rowranges <- GRanges("1", IRanges(1 + 0:3 * 10, width=10))
se <- SummarizedExperiment(list(counts=counts), colData=coldata, rowRanges=rowranges)
mcols(colData(se)) <- DataFrame(info="x is a factor")
mcols(se)$id <- 1:4
mcols(mcols(se)) <- DataFrame(info="the gene id")
dds <- DESeqDataSet(se, ~ x)
mcols(colData(dds))
mcols(mcols(dds))

