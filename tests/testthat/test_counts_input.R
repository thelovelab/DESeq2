# count matrix input
cnts <- matrix(rnbinom(40,mu=100,size=2),ncol=4)
mode(cnts) <- "integer"
coldata <- data.frame(cond=factor(c("A","A","B","B")))
dds <- DESeqDataSetFromMatrix(cnts, coldata, ~cond)

# error when assay colnames not equal to colData rownames
colnames(cnts) <- letters[1:4]
expect_error(DESeqDataSetFromMatrix(cnts, coldata, ~cond))

# tidy data frame input
gene.names <- paste0("gene",1:10)
tidy.counts <- cbind(gene.names, as.data.frame(cnts))
rownames(coldata) <- letters[1:4]
dds <- DESeqDataSetFromMatrix(tidy.counts, coldata, ~cond, tidy=TRUE)
expect_true(all(rownames(dds) == gene.names))
expect_true(all(counts(dds) == cnts))

