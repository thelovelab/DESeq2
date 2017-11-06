context("design matrix")
test_that("design can be a matrix", {

  m <- matrix(rpois(1000,100),ncol=12)
  coldata <- data.frame(condition=factor(rep(1:2,each=6)),
                        batch=factor(rep(c(1,2,1,2),each=3)))
  dm <- model.matrix(~condition, coldata)
  dm2 <- model.matrix(~batch + condition, coldata)
  dds <- DESeqDataSetFromMatrix(m, coldata, dm)
  dds <- DESeq(dds, fitType="mean")
  resultsNames(dds)
  dds <- DESeq(dds, full=dm2, fitType="mean")
  resultsNames(dds)
  
})
