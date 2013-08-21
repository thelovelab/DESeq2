# test the optimization of the logarithm of dispersion (alpha)
# parameter with Cox-Reid adjustment and prior distribution.
# also test the derivatives of the log posterior w.r.t. log alpha
test_dispersionFitting <- function() {
  m <- 10
  set.seed(1)
  y <- rpois(m,20)
  sf <- rep(1,m)
  x <- cbind(rep(1,m),rep(0:1,each=m/2))
  colnames(x) <- c("Intercept","condition")
  
  lambda <- 2
  alpha <- .5

  # make a DESeqDataSet but don't use the design formula
  # instead we supply a model matrix below
  dds <- DESeqDataSetFromMatrix(matrix(y,nrow=1),
                                colData=DataFrame(condition=x[,2]),
                                design= ~ condition)
  sizeFactors(dds) <- sf
  dispersions(dds) <- alpha
  mcols(dds)$baseMean <- mean(y)

  # for testing we convert beta to the naturual log scale:
  # convert lambda from log to log2 scale by multiplying by log(2)^2
  # then convert beta back from log2 to log scale by multiplying by log(2)
  betaDESeq <- log(2)*DESeq2:::fitNbinomGLMs(dds, lambda=c(0,lambda*log(2)^2),
                                             modelMatrix=x)$betaMatrix

  log_alpha_prior_mean <- .5
  log_alpha_prior_sigmasq <- 1
  mu.hat <- as.numeric(exp(x %*% t(betaDESeq)))
  
  dispRes <- DESeq2:::fitDisp(ySEXP = matrix(y,nrow=1), xSEXP = x,
                              mu_hatSEXP = matrix(mu.hat,nrow=1), log_alphaSEXP = 0,
                              log_alpha_prior_meanSEXP = log_alpha_prior_mean,
                              log_alpha_prior_sigmasqSEXP = log_alpha_prior_sigmasq,
                              min_log_alphaSEXP = log(1e-8), kappa_0SEXP = 1,
                              tolSEXP = 1e-16, maxitSEXP = 100, use_priorSEXP = TRUE)
  
  # maximum a posteriori (MAP) estimate from DESeq
  dispDESeq <- dispRes$log_alpha

  # MAP estimate using optim
  logPost <- function(log.alpha) {
    alpha <- exp(log.alpha)
    w <- diag(1/(1/mu.hat^2 * ( mu.hat + alpha * mu.hat^2 )))
    logLike <- sum(dnbinom(y, mu=mu.hat, size=1/alpha, log=TRUE))
    coxReid <- -.5*(log(det(t(x) %*% w %*% x)))
    prior <- dnorm(log.alpha, log_alpha_prior_mean, sqrt(log_alpha_prior_sigmasq), log=TRUE)
    (logLike + coxReid + prior)
  }

  dispOptim <- optim(0, function(p) -1*logPost(p), control=list(reltol=1e-16),
                     method="Brent", lower=-10, upper=10)$par

  checkEqualsNumeric(dispDESeq, dispOptim, tolerance=1e-6)
    
  # check derivatives:
  
  # from Ted Harding https://stat.ethz.ch/pipermail/r-help/2007-September/140013.html
  num.deriv <- function(f,x,h=0.001) (f(x + h/2) - f(x-h/2))/h
  num.2nd.deriv <- function(f,x,h=0.001) (f(x + h) - 2*f(x) + f(x - h))/h^2
  
  # first derivative of log posterior w.r.t log alpha at start
  dispDerivDESeq <- dispRes$initial_dlp
  dispDerivNum <- num.deriv(logPost,0)

  checkEqualsNumeric(dispDerivDESeq, dispDerivNum, tolerance=1e-6)
  
  # second derivative at finish
  dispD2DESeq <- dispRes$last_d2lp
  dispD2Num <- num.2nd.deriv(logPost, dispRes$log_alpha)

  checkEqualsNumeric(dispD2DESeq, dispD2Num, tolerance=1e-6)
}


test_alternativeDispersions <- function() {
  dds <- makeExampleDESeqDataSet()
  dds <- estimateSizeFactors(dds)
  
  ddsLocal <- estimateDispersions(dds, fitType="local")
  
  ddsMean <- estimateDispersions(dds, fitType="mean")
  
  ddsMed <- estimateDispersionsGeneEst(dds)
  useForMedian <- mcols(ddsMed)$dispGeneEst > 1e-7
  medianDisp <- median(mcols(ddsMed)$dispGeneEst[useForMedian],na.rm=TRUE)
  mcols(ddsMed)$dispFit <- medianDisp
  ddsMed <- estimateDispersionsMAP(ddsMed)
}
