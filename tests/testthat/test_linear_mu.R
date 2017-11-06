context("linear_mu")
test_that("the use of linear model for fitting mu works as expected", {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100, m=4, interceptMean=10, interceptSD=3,
                                 dispMeanRel=function(x) 0.5, sizeFactors=c(.5,1,1,2))
  dds <- estimateSizeFactors(dds)
  dds1 <- estimateDispersionsGeneEst(dds, linearMu=FALSE)
  dds2 <- estimateDispersionsGeneEst(dds, linearMu=TRUE)
  mu1 <- assays(dds1)[["mu"]]
  mu2 <- assays(dds2)[["mu"]]
  #par(mfrow=c(2,2),mar=c(3,3,1,1))
  #for (i in 1:4) {
    #plot(mu1[,i], mu2[,i], xlab="", ylab="", log="xy")
    #abline(0,1)
  #}
  cors <- diag(cor(mu1, mu2, use="complete"))
  expect_true(all(cors > 1 - 1e-6))
  #
  dds2 <- estimateDispersionsFit(dds2, fitType="mean")
  dds2 <- estimateDispersionsMAP(dds2)
  dds2 <- nbinomWaldTest(dds2)
  res <- results(dds2)
})
