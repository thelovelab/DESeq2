context("lfcShrink")
test_that("LFC shrinkage works", {
  
  # testing out various methods for LFC shrinkage
  set.seed(1)
  dds <- makeExampleDESeqDataSet(betaSD=1,n=1000,m=20)
  expect_error(lfcShrink(dds, coef=2), "first run")
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

  # LFC threshold for "normal" and "apeglm"
  res0 <- results(dds, name="condition_B_vs_A", lfcThreshold=1)
  res.n <- lfcShrink(dds=dds, coef=2, type="normal", lfcThreshold=1)
  res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm", lfcThreshold=1)
  
  #plotMA(res0, ylim=c(-4,4), cex=1); abline(h=c(-1,1),col="dodgerblue")
  #plotMA(res.n, ylim=c(-4,4), cex=1); abline(h=c(-1,1),col="dodgerblue")
  #plotMA(res.ape, ylim=c(-4,4), cex=1); abline(h=c(-1,1),col="dodgerblue")
  
  # s-value returned
  res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm", svalue=TRUE)
  expect_true("svalue" %in% names(res.ape))
  res.ash <- lfcShrink(dds=dds, res=res, type="ashr", svalue=TRUE)
  expect_true("svalue" %in% names(res.ash))

  # plotMA works with s-values
  plotMA(res.ape, cex=1)
  plotMA(res.ash, cex=1)

  # summary works with s-values
  summary(res.ape)
  summary(res.ash)
  
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
  res.normal <- lfcShrink(dds=dds, coef=2, res=res, type="normal")
  res.ape <- lfcShrink(dds=dds, coef=2, res=res, type="apeglm")
  
})  
