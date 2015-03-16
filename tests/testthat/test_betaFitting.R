# test for equivalence of DESeq2 estimates with those
# found using IRLS code and using optim
m <- 10
set.seed(1)
y <- rpois(m,20)
sf <- rep(1,m)
condition <- factor(rep(0:1,each=m/2))
x <- cbind(rep(1,m),rep(0:1,each=m/2))
lambda <- 2
alpha <- .5

dds <- DESeqDataSetFromMatrix(matrix(y,nrow=1),
                              colData=DataFrame(condition),
                              design= ~ condition)
sizeFactors(dds) <- sf
dispersions(dds) <- alpha
mcols(dds)$baseMean <- mean(y)

# for testing we convert beta to the naturual log scale:
# convert lambda from log to log2 scale by multiplying by log(2)^2
# then convert beta back from log2 to log scale by multiplying by log(2)
betaDESeq <- log(2)*DESeq2:::fitNbinomGLMs(dds, lambda=c(0,lambda*log(2)^2))$betaMatrix

# the IRLS algorithm
betaIRLS <- c(1,1)
for (t in 1:100) {
  mu.hat <- as.vector(sf * exp(x %*% betaIRLS))
  w <- diag(1/(1/mu.hat^2 * ( mu.hat + alpha * mu.hat^2 )))
  z <- log(mu.hat/sf) + (y - mu.hat)/mu.hat
  ridge <- diag(c(0,lambda))
  betaIRLS <- as.vector(solve(t(x) %*% w %*% x + ridge) %*% t(x) %*% w %*% z)
}

# using optim
objectiveFn <- function(p) {
  mu <- exp(x %*% p)
  logLike <- sum(dnbinom(y, mu=mu, size=1/alpha, log=TRUE))
  prior <- dnorm(p[2], 0, sqrt(1/lambda),log=TRUE)
  -1 * (logLike + prior)
}
betaOptim <- optim(c(.1,.1), objectiveFn, control=list(reltol=1e-16))$par

expect_equal(as.numeric(betaDESeq), betaIRLS, tolerance=1e-6)
expect_equal(as.numeric(betaDESeq), betaOptim, tolerance=1e-6)
