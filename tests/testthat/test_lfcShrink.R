context("lfcShrink")
test_that("LFC shrinkage works", {
  dds <- makeExampleDESeqDataSet(betaSD=1)
  dds <- estimateSizeFactors(dds)
  expect_error(lfcShrink(dds, 2, 1))
  dds <- estimateDispersions(dds)
  dds <- DESeq(dds)
  res <- results(dds)
  res.shr <- lfcShrink(dds=dds, coef=2, res=res)
  #plotMA(res.shr)
  res.shr <- lfcShrink(dds=dds,
                       contrast=c("condition","B","A"),
                       res=res)
  #plotMA(res.shr)

  # testing out various methods for LFC shrinkage
  set.seed(1)
  dds <- makeExampleDESeqDataSet(betaSD=1,n=1000,m=10)
  dds <- DESeq(dds)
  res <- results(dds, name="condition_B_vs_A")

  # dds and res must match
  expect_error(lfcShrink(dds=dds, coef=2, res=res[1:500,], type="normal"), "rownames")
  expect_error(lfcShrink(dds=dds, coef=2, res=res[1:500,], type="apeglm"), "rownames")  

  # try out various types and ways of specifying coefs
  res.n <- lfcShrink(dds=dds, coef="condition_B_vs_A", res=res, type="normal")
  res.n <- lfcShrink(dds=dds, coef=2, res=res, type="normal")
  res.n <- lfcShrink(dds=dds, coef=2, type="normal")
  res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm")
  res.ash <- lfcShrink(dds=dds, res=res, type="ashr")

  # prior info
  ## str(priorInfo(res.n))
  ## str(priorInfo(res.ape))
  ## str(priorInfo(res.ash))

  # plot against true
  ## par(mfrow=c(1,3))
  ## plot(mcols(dds)$trueBeta, res.n$log2FoldChange); abline(0,1,col="red")
  ## plot(mcols(dds)$trueBeta, res.ape$log2FoldChange); abline(0,1,col="red")
  ## plot(mcols(dds)$trueBeta, res.ash$log2FoldChange); abline(0,1,col="red")
  
  # s-value returned
  res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm", svalue=TRUE)
  expect_true("svalue" %in% names(res.ape))
  res.ash <- lfcShrink(dds=dds, res=res, type="ashr", svalue=TRUE)
  expect_true("svalue" %in% names(res.ash))

  # TODO add tests of new plotMA() with svalue
  
  # list returned
  res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm", returnList=TRUE)
  names(res.ape)
  res.ash <- lfcShrink(dds=dds, res=res, type="ashr", returnList=TRUE)
  names(res.ash)

  # test wrong coef specified
  resInt <- results(dds, name="Intercept")
  expect_error(lfcShrink(dds=dds, coef=2, res=resInt, type="apeglm"))

  # test supplied model.matrix
  full <- model.matrix(~condition, colData(dds))
  dds <- DESeq(dds, full=full)
  res <- results(dds)
  res.ape <- lfcShrink(dds=dds, coef=2, res=res, type="apeglm")
  
})  
