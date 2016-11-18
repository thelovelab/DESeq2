context("lfcShrink")
test_that("LFC shrinkage works", {
  dds <- makeExampleDESeqDataSet(betaSD=1)
  dds <- estimateSizeFactors(dds)
  expect_error(lfcShrink(dds, 2, 1))
  dds <- estimateDispersions(dds)
  lfc <- lfcShrink(dds=dds, coef=2)
  dds <- DESeq(dds, betaPrior=FALSE)
  res <- results(dds)
  res.shr <- lfcShrink(dds=dds, coef=2, res=res)
  plotMA(res.shr)
  res.shr <- lfcShrink(dds=dds,
                       contrast=c("condition","B","A"),
                       res=res)
  plotMA(res.shr)
})  
