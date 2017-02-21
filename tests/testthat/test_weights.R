context("weights")
test_that("weights work", {
  dds <- makeExampleDESeqDataSet()
  dds <- DESeq(dds)
  dds2 <- dds
  assays(dds2)[["weights"]] <- matrix(1, nrow=1000, ncol=12)
  assays(dds2)[["weights"]][1,1] <- 0
  dds2 <- nbinomWaldTest(dds2)
  dds3 <- dds[,-1]
  dds3 <- nbinomWaldTest(dds3)
  rbind(results(dds)[1,],
        results(dds2)[1,],
        results(dds3)[1,])

  # need to test intercept code
  # need to test optim code

  # need to implement dispersion weights
  
})
