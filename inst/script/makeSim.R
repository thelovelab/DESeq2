makeSim <- function(n, m, x, beta, meanDispPairs, sf=rep(1,m)) {
  idx <- sample(nrow(meanDispPairs), n, replace=TRUE)
  mu0 <- meanDispPairs[idx,1]
  disp <- meanDispPairs[idx,2]
  betafull <- cbind(log2(mu0), beta)
  mu <- 2^(betafull %*% t(x))
  muMat <- matrix(rep(mu, times=m) * rep(sf, each=n), ncol=m)
  list(mat = matrix(rnbinom(n*m, mu=muMat, size=1/disp), ncol=m),
       disp = disp,
       mu0 = mu0)
}
