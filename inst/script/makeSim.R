makeSim <- function(n, m, x, beta, meanDispPairs, sf=rep(1,m)) {
  idx <- sample(nrow(meanDispPairs), n, replace=TRUE)
  mu0 <- meanDispPairs[idx,1]
  disp <- meanDispPairs[idx,2]
  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  swap <- sort(sample(n, n/2, replace=FALSE))
  mu[swap,] <- mu[swap,c((m/2+1):m,1:(m/2))]
  beta[swap] <- -1 * beta[swap]
  muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)
  list(mat = matrix(rnbinom(n*m, mu=muMat, size=1/disp), ncol=m),
       beta = beta,
       disp = disp,
       mu0 = mu0)
}
