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
  dds <- DESeq(dds, full=dm2, fitType="mean")
  resultsNames(dds)

  res <- results(dds)
  res <- results(dds, contrast=list("condition2","batch2"))
  res <- results(dds, contrast=c(0,-1,1))
  expect_error(res <- results(dds, contrast=c("condition","2","1")), "only list- and numeric-type")

  # TODO make this work with type="normal"
  #res <- lfcShrink(dds, coef="condition2", type="normal")
  res <- lfcShrink(dds, coef="condition2", type="apeglm")
  res <- lfcShrink(dds, coef="condition2", type="ashr")
  
  # test replace with matrix
  design(dds) <- dm
  
})
