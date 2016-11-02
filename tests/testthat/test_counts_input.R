context("counts_input")
test_that("counts can be supplied as input (tidy or not)", {
  # count matrix input
  cnts <- matrix(rnbinom(40,mu=100,size=2),ncol=4)
  mode(cnts) <- "integer"
  coldata <- data.frame(cond=factor(c("A","A","B","B")))
  dds <- DESeqDataSetFromMatrix(cnts, coldata, ~cond)

  # tidy data frame input
  gene.names <- paste0("gene",1:10)
  rownames(coldata) <- colnames(cnts) <- letters[1:4]
  tidy.counts <- cbind(gene.names, as.data.frame(cnts))
  dds <- DESeqDataSetFromMatrix(tidy.counts, coldata, ~cond, tidy=TRUE)
  expect_true(all(rownames(dds) == gene.names))
  expect_true(all(counts(dds) == cnts))
})
