context("no replicates")
test_that("no replicates gives deprecation and usage warning", {
  
  dds <- makeExampleDESeqDataSet(n=100, m=2)
  expect_warning({ dds <- DESeq(dds) })
  res <- results(dds)
  
})
