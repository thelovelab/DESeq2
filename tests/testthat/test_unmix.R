context("unmix")
test_that("unmixing samples works", {

  set.seed(1)
  n <- 2000
  a <- runif(n)
  b <- runif(n)
  c <- runif(n)
  counts <- matrix(nrow=n, ncol=8)
  disp <- 0.01

  counts[,1] <- rnbinom(n, mu=1e4 * a, size=1/disp)
  counts[,2] <- rnbinom(n, mu=1e4 * b, size=1/disp)
  counts[,3] <- rnbinom(n, mu=1e4 * c, size=1/disp)
  counts[,4] <- rnbinom(n, mu=1e4 * (.75*a + .25*b), size=1/disp)
  counts[,5] <- rnbinom(n, mu=1e4 * (.5*a + .5*b), size=1/disp)
  counts[,6] <- rnbinom(n, mu=1e4 * (.25*a + .75*b), size=1/disp)
  counts[,7] <- rnbinom(n, mu=1e4 * (.33*a + .33*b + .33*c), size=1/disp)
  counts[,8] <- rnbinom(n, mu=1e4 * (.25*a + .25*b + .5*c), size=1/disp)
  coldata <- data.frame(a=c(1,0,0,.75,.5,.25,.33,.25),
                        b=c(0,1,0,.25,.5,.75,.33,.25),
                        c=c(0,0,1,  0,  0, 0,.33,.5))

  dds <- DESeqDataSetFromMatrix(counts, coldata, ~1)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType="mean")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  #library(ggplot2)
  #plotPCA(vsd, intgroup=c("a","b","c")) + geom_point(size=5)

  pure <- matrix(rnbinom(3*n,mu=1e4*c(a,b,c),size=1/disp),ncol=3)
  colnames(pure) <- c("a","b","c")

  colnames(counts) <- paste0("sample",seq_len(ncol(counts)))
  
  alpha <- attr(dispersionFunction(dds),"mean")

  mix <- unmix(counts, pure=pure, alpha=alpha, quiet=TRUE)

  max(abs(dds$a - mix[,1]))
  max(abs(dds$b - mix[,2]))
  max(abs(dds$c - mix[,3]))
  
  expect_lt(max(abs(dds$a - mix[,1])), .01)
  expect_lt(max(abs(dds$b - mix[,2])), .01)
  expect_lt(max(abs(dds$c - mix[,3])), .01)

  # test the shifted log (designed for TPMs)
  mix2 <- unmix(counts, pure=pure, shift=0.5, quiet=TRUE)

  # test expanded output
  mix3 <- unmix(counts, pure=pure, alpha=alpha, format="list", quiet=TRUE)

  # test warning
  pure <- cbind(pure[,1:3], d=(pure[,3] + rpois(n, 10)))
  expect_warning(unmix(counts, pure=pure, alpha=alpha, quiet=TRUE), "highly correlated")
  
})
