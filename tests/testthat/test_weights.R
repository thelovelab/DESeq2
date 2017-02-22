context("weights")
test_that("weights work", {
  set.seed(1)
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

  mcols(dds2)[1,"deviance"]
  mcols(dds3)[1,"deviance"]
  
  nf <- matrix(sizeFactors(dds),nrow=nrow(dds),ncol=ncol(dds),byrow=TRUE)
  
  o <- fitNbinomGLMsOptim(object=dds,
                          modelMatrix=model.matrix(design(dds), colData(dds)),
                          lambda=rep(1e-6, 2),
                          rowsForOptim=1,
                          rowStable=TRUE,
                          normalizationFactors=nf,
                          alpha_hat=dispersions(dds),
                          weights=w,
                          useWeights=TRUE,
                          betaMatrix=matrix(0,nrow=nrow(dds),ncol=2),
                          betaSE=matrix(0,nrow=nrow(dds),ncol=2),
                          betaConv=rep(FALSE,nrow(dds)),
                          beta_mat=matrix(0,nrow=nrow(dds),ncol=2),
                          mu=matrix(0,nrow=nrow(dds),ncol=ncol(dds)),
                          logLike=rep(0,nrow(dds)))

  results(dds3)[1,2:3]
  c(o$betaMatrix[1,2], o$betaSE[1,2])
  
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

  mcols(dds2)[1,"deviance"]
  mcols(dds3)[1,"deviance"]
  
  # need to implement dispersion weights

  # some code for testing whether still full rank
  x <- model.matrix(~ condition, colData(dds))
  w <- matrix(runif(20000*12), nrow=20000, ncol=12)
  test <- logical(nrow(w))
  m <- ncol(x)
  system.time({
    for (i in 1:nrow(w)) {
      test[i] <- qr(w[i,] * x)$rank == m
    }
  })
  stopifnot(all(test))
  
})
