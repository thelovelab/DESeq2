context("weights")
test_that("weights work", {
  dds <- makeExampleDESeqDataSet(n=100)
  dds <- DESeq(dds, quiet=TRUE)
  dds2 <- dds
  w <- matrix(1, nrow=nrow(dds), ncol=12)
  w[1,1] <- 0
  assays(dds2)[["weights"]] <- w
  dds2 <- nbinomWaldTest(dds2)
  dds3 <- dds[,-1]
  dds3 <- nbinomWaldTest(dds3)
  rbind(results(dds)[1,],
        results(dds2)[1,],
        results(dds3)[1,])

  o <- fitNbinomGLMsOptim(object=dds,
                          modelMatrix=model.matrix(design(dds), colData(dds)),
                          lambda=rep(1e-6, 2),
                          rowsForOptim=1,
                          rowStable=TRUE,
                          normalizationFactors=matrix(sizeFactors(dds),
                                                      nrow=nrow(dds),
                                                      ncol=ncol(dds),
                                                      byrow=TRUE),
                          alpha_hat=dispersions(dds),
                          weights=w,
                          useWeights=TRUE,
                          betaMatrix=matrix(0,nrow=nrow(dds),ncol=2),
                          betaSE=matrix(0,nrow=nrow(dds),ncol=2),
                          betaConv=rep(FALSE,nrow(dds)),
                          beta_mat=matrix(0,nrow=nrow(dds),ncol=2),
                          mu=matrix(0,nrow=nrow(dds),ncol=ncol(dds)),
                          logLike=rep(0,nrow(dds)))

  results(dds3)[1,2]
  o$betaMatrix[1,2]
  
  design(dds) <- ~1
  suppressWarnings({ dds <- DESeq(dds, quiet=TRUE) })
  dds2 <- dds
  assays(dds2)[["weights"]] <- w
  dds2 <- nbinomWaldTest(dds2)
  dds3 <- dds[,-1]
  dds3 <- nbinomWaldTest(dds3)
  rbind(results(dds)[1,],
        results(dds2)[1,],
        results(dds3)[1,])

  # need to implement dispersion weights
  
})
