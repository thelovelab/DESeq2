context("lfcShrink")
test_that("LFC shrinkage works", {
  dds <- makeExampleDESeqDataSet(betaSD=1)
  dds <- estimateSizeFactors(dds)
  expect_error(lfcShrink(dds, 2, 1))
  dds <- estimateDispersions(dds)
  lfc <- lfcShrink(dds=dds, coef=2)
  dds <- DESeq(dds)
  res <- results(dds)
  res.shr <- lfcShrink(dds=dds, coef=2, res=res)
  plotMA(res.shr)
  res.shr <- lfcShrink(dds=dds,
                       contrast=c("condition","B","A"),
                       res=res)
  plotMA(res.shr)

  # testing out various methods for LFC shrinkage
  set.seed(1)
  dds <- makeExampleDESeqDataSet(betaSD=1,n=1000,m=10)
  dds <- dds[rowSums(counts(dds)) > 0,]
  dds <- DESeq(dds)
  res <- results(dds, name="condition_B_vs_A")
  res.n <- lfcShrink(dds=dds, coef=2, res=res, type="normal")
  res.ape <- lfcShrink(dds=dds, coef=2, res=res, type="apeglm")
  res.ash <- lfcShrink(dds=dds, res=res, type="ashr")

  str(priorInfo(res.n))
  str(priorInfo(res.ape))
  str(priorInfo(res.ash))

  par(mfrow=c(1,3))
  plot(mcols(dds)$trueBeta, res.n$log2FoldChange); abline(0,1,col="red")
  plot(mcols(dds)$trueBeta, res.ape$log2FoldChange); abline(0,1,col="red")
  plot(mcols(dds)$trueBeta, res.ash$log2FoldChange); abline(0,1,col="red")
  
})  
