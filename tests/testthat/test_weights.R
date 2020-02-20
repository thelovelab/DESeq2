context("weights")
test_that("weights work", {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=10)

  # check that weight to 0 is like a removed samples
  dds <- DESeq(dds, quiet=TRUE)
  dds2 <- dds
  w <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
  w[1,1] <- 0
  assays(dds2, withDimnames=FALSE)[["weights"]] <- w
  dds2 <- nbinomWaldTest(dds2)
  dds3 <- dds[,-1]
  dds3 <- nbinomWaldTest(dds3)

  # in terms of LFC, SE and deviance
  expect_equal(results(dds2)$log2FoldChange[1], results(dds3)$log2FoldChange[1])
  expect_equal(results(dds2)$lfcSE[1], results(dds3)$lfcSE[1])
  expect_equal(mcols(dds2)[1,"deviance"],mcols(dds3)[1,"deviance"])

  # check weights working in the optim code
  
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

  expect_equal(results(dds3)$log2FoldChange[1], o$betaMatrix[1,2], tolerance=1e-4)

  # can use weights with betaPrior=TRUE
  
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=10)
  w <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
  w[1,1] <- 0
  assays(dds, withDimnames=FALSE)[["weights"]] <- w
  dds <- DESeq(dds, betaPrior=TRUE, quiet=TRUE)

  # check weights working for intercept only
  
  design(dds) <- ~1
  suppressWarnings({ dds <- DESeq(dds, quiet=TRUE) })
  dds2 <- dds
  assays(dds2, withDimnames=FALSE)[["weights"]] <- w
  dds2 <- nbinomWaldTest(dds2)
  dds3 <- dds[,-1]
  dds3 <- nbinomWaldTest(dds3)

  expect_equal(results(dds2)$log2FoldChange[1], results(dds3)$log2FoldChange[1])
  expect_equal(results(dds2)$lfcSE[1], results(dds3)$lfcSE[1])
  expect_equal(mcols(dds2)[1,"deviance"],mcols(dds3)[1,"deviance"])

  # check that weights downweight outlier in dispersion estimation
  
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=10)
  counts(dds)[1,1] <- 100L
  sizeFactors(dds) <- rep(1,ncol(dds))
  dds <- estimateDispersions(dds)
  dds2 <- dds
  w <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
  w[1,1] <- 0
  assays(dds2, withDimnames=FALSE)[["weights"]] <- w
  dds2 <- estimateDispersions(dds2)
  dds3 <- dds[,-1]
  dds3 <- estimateDispersions(dds3)
  
  expect_equal(mcols(dds2)[1,"dispGeneEst"],mcols(dds3)[1,"dispGeneEst"],tolerance=1e-3)
  # MAP estimates won't be equal because of different dispersion prior widths...
  expect_true(mcols(dds)[1,"dispMAP"] > mcols(dds2)[1,"dispMAP"])
  
})

test_that("weights failing check gives warning, passes them through", {

  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=10)
  w <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
  w[1,1:6] <- 0
  assays(dds, withDimnames=FALSE)[["weights"]] <- w
  expect_warning(dds <- DESeq(dds))
  expect_true(mcols(dds)$allZero[1])
  expect_true(mcols(dds)$weightsFail[1])
  res <- results(dds)
  
})

test_that("weights with and without CR term included", {

  set.seed(1); alpha <- .25
  dmr <- function(x) alpha
  dds <- makeExampleDESeqDataSet(n=50, m=100, betaSD=1, interceptMean=10, interceptSD=.5, dispMeanRel=dmr)
  dds$group <- factor(rep(1:50,2)); design(dds) <- ~0 + group + condition
  w <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
  o <- 35
  w[,c(1:o, 50 + 1:o)] <- 1e-6
  assays(dds, withDimnames=FALSE)[["weights"]] <- w
  counts(dds)[,c(1:o, 50 + 1:o)] <- 1L
  sizeFactors(dds) <- 1
  dds <- estimateDispersions(dds, fitType="mean")
  dds2 <- estimateDispersions(dds, fitType="mean", useCR=FALSE)
  ## dds3 <- estimateDispersions(dds, fitType="mean", weightThreshold=0)
  ## par(mfcol=c(3,2), mar=c(4.5,4.5,1,1), cex=1.1)
  ## for (col in c("dispGeneEst","dispersion")) {
  ##   plot(mcols(dds)$trueBeta, mcols(dds)[[col]], ylim=c(0,2*alpha));abline(h=alpha,col="blue")
  ##   plot(mcols(dds2)$trueBeta, mcols(dds2)[[col]], ylim=c(0,2*alpha));abline(h=alpha,col="blue")
  ##   plot(mcols(dds3)$trueBeta, mcols(dds3)[[col]], ylim=c(0,10));abline(h=alpha,col="blue")
  ## }
  
})
