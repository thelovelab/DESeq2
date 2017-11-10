context("design matrix")
test_that("design can be a matrix", {

  m <- matrix(rpois(12*100,100),ncol=12)
  coldata <- data.frame(condition=factor(rep(1:2,each=6)),
                        batch=factor(rep(c(1,2,1,2),each=3)))
  dm <- model.matrix(~condition, coldata)
  dm2 <- model.matrix(~batch + condition, coldata)
  dds <- DESeqDataSetFromMatrix(m, coldata, dm)
  dds <- DESeq(dds, fitType="mean")
  resultsNames(dds)
  
  # specifying 'full' overrides...
  dds2 <- DESeq(dds, full=dm2, fitType="mean")
  resultsNames(dds2)

  dds <- DESeqDataSetFromMatrix(m, coldata, dm2)
  dds <- DESeq(dds, fitType="mean")
  
  res <- results(dds)
  res <- results(dds, contrast=list("condition2","batch2"))
  res <- results(dds, contrast=c(0,-1,1))
  expect_error(res <- results(dds, contrast=c("condition","2","1")), "only list- and numeric-type")

  res <- lfcShrink(dds, coef="condition2", type="normal")
  res <- lfcShrink(dds, coef="condition2", type="apeglm")
  res <- lfcShrink(dds, coef="condition2", type="ashr")
  
  # test replace with matrix
  design(dds) <- dm
  
})
