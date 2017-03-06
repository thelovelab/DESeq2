context("size_factor")
test_that("size factor works", {
  # size factor tests
  m <- matrix(1:16, ncol=4)
  expect_error(estimateSizeFactorsForMatrix(m, geoMeans=1:5))
  expect_error(estimateSizeFactorsForMatrix(m, geoMeans=rep(0,4)))
  expect_error(estimateSizeFactorsForMatrix(m, controlGenes="foo"))
  estimateSizeFactorsForMatrix(m, geoMeans=1:4)
  estimateSizeFactorsForMatrix(m, controlGenes=1:2)
  
  # iterate method
  set.seed(1)
  true.sf <- 2^(rep(c(-2,-1,0,0,1,2),each=2))
  dmr <- function(x) 0.01
  dds <- makeExampleDESeqDataSet(sizeFactors=true.sf, n=100, dispMeanRel=dmr)
  cts <- counts(dds)
  idx <- cbind(seq_len(nrow(cts)), sample(ncol(dds), nrow(cts), replace=TRUE))
  cts[idx] <- 0L
  cts[1,1] <- 1000000L # an outlier
  counts(dds) <- cts

  dds <- estimateSizeFactors(dds, type="iterate")
  sf <- sizeFactors(dds)
  plot(true.sf, sf);abline(0,1)
  coefs <- coef(lm(sf ~ true.sf))
  expect_true(abs(coefs[1]) < .1)
  expect_true(abs(coefs[2] - 1) < .1)

  dds <- estimateSizeFactors(dds, type="poscounts")
  sf <- sizeFactors(dds)
  plot(true.sf, sf);abline(0,1)
  coefs <- coef(lm(sf ~ true.sf))
  expect_true(abs(coefs[1]) < .1)
  expect_true(abs(coefs[2] - 1) < .1)

})
