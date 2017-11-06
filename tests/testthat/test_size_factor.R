context("size_factor")
test_that("size factor works", {
  
  # size factor error checking
  m <- matrix(1:16, ncol=4)
  expect_error(estimateSizeFactorsForMatrix(m, geoMeans=1:5))
  expect_error(estimateSizeFactorsForMatrix(m, geoMeans=rep(0,4)))
  expect_error(estimateSizeFactorsForMatrix(m, controlGenes="foo"))
  estimateSizeFactorsForMatrix(m, geoMeans=1:4)
  estimateSizeFactorsForMatrix(m, controlGenes=1:2)

  # norm matrix works
  nm <- m / exp(rowMeans(log(m))) # divide out the geometric mean
  true.sf <- c(2,1,1,.5)
  counts <- sweep(2*m, 2, true.sf, "*")
  dds <- DESeqDataSetFromMatrix(counts, data.frame(x=1:4), ~1)
  dds <- estimateSizeFactors(dds, normMatrix=nm)
  expect_equal((normalizationFactors(dds)/nm)[1,], true.sf)
  
  # make some counts with zeros
  set.seed(1)
  true.sf <- 2^(rep(c(-2,-1,0,0,1,2),each=2))
  dmr <- function(x) 0.01
  dds <- makeExampleDESeqDataSet(sizeFactors=true.sf, n=100, dispMeanRel=dmr)
  cts <- counts(dds)
  idx <- cbind(seq_len(nrow(cts)), sample(ncol(dds), nrow(cts), replace=TRUE))
  cts[idx] <- 0L
  cts[1,1] <- 1000000L # an outlier
  counts(dds) <- cts

  # positive counts method
  dds <- estimateSizeFactors(dds, type="poscounts")
  sf <- sizeFactors(dds)
  #plot(true.sf, sf);abline(0,1)
  coefs <- coef(lm(sf ~ true.sf))
  expect_true(abs(coefs[1]) < .1)
  expect_true(abs(coefs[2] - 1) < .1)
  
  # iterate method
  dds <- estimateSizeFactors(dds, type="iterate")
  sf <- sizeFactors(dds)
  #plot(true.sf, sf);abline(0,1)
  coefs <- coef(lm(sf ~ true.sf))
  expect_true(abs(coefs[1]) < .1)
  expect_true(abs(coefs[2] - 1) < .1)

})
