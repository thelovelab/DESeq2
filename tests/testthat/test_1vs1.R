context("1vs1")
test_that("1 vs 1 gets warning", {
  dds <- makeExampleDESeqDataSet(n=100, m=2)
  expect_warning({ dds <- DESeq(dds)})
  res <- results(dds)
})
