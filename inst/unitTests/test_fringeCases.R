test_oneRow <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=1)
  sizeFactors(dds) <- rep(1,ncol(dds))
  dispersions(dds) <- .5
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  dds <- nbinomLRT(dds, reduced=~1)
  res <- results(dds)
}

test_onlyIntercept <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100)
  design(dds) <- ~ 1
  dds <- DESeq(dds)
  res <- results(dds)
}
