context("interactions")
test_that("interactions throw error", {
  dds <- makeExampleDESeqDataSet(n=100,m=8)
  colData(dds)$group <- factor(rep(c("X","Y"),times=ncol(dds)/2))
  design(dds) <- ~ condition + group + condition:group
  dds <- DESeq(dds)
  expect_equal(resultsNames(dds)[4], "conditionB.groupY")
  # interactions error
  expect_error(DESeq(dds, betaPrior=TRUE), "designs with interactions")

  # also lfcShrink
  res <- results(dds, name="conditionB.groupY")
  expect_error(res <- lfcShrink(dds, coef=4, res=res), "not implemented")

  res <- results(dds, contrast=c("condition","B","A"))
  expect_error(res <- lfcShrink(dds, contrast=c("condition","B","A"), res=res),
               "not implemented")  
})
