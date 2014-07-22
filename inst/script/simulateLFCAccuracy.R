source("makeSim.R")
source("runScripts.R")
load("meanDispPairs.RData")
library("DESeq2")
library("edgeR")
library("Biobase")
algos <- list("DESeq2"=runDESeq2,"edgeR"=runEdgeR)
namesAlgos <- names(algos)

set.seed(1)
nreps <- 10
n <- 1000
ms <- c(4,6,10,16,20)
types <- c("bell","slab bell","slab spike","spike spike")
methods <- c("DESeq2","edgeR predFC","edgeR predFC10")
res <- do.call(rbind, lapply(ms, function(m) {
  do.call(rbind, lapply(types, function(type) {
    do.call(rbind, lapply(seq_len(nreps), function(i) {
      beta <- if (type == "bell") {
        rnorm(n)
      } else if (type == "slab bell") {
        c(rnorm(n * 8/10), runif(n * 2/10, -4, 4))
      } else if (type == "slab spike") {
        beta <- c(rep(0, n * 8/10), runif(n * 2/10, -4, 4))
      } else if (type == "spike spike") {
        beta <- c(rep(0, n * 8/10), sample(c(-2, 2), n * 2/10, TRUE))
      }
      condition <- factor(rep(c("A","B"), each = m/2))
      x <- model.matrix(~ condition)
      mat <- makeSim(n,m,x,beta,meanDispPairs)$mat      
      e <- ExpressionSet(mat, AnnotatedDataFrame(data.frame(condition)))

      res0 <- lapply(algos, function(f) f(e))
      resdf <- do.call(rbind, lapply(methods, function(method) {
        if (method == "edgeR predFC") {
          predbeta <- res0[["edgeR"]]$predbeta
        } else if (method == "edgeR predFC10") {
          predbeta <- res0[["edgeR"]]$predbeta10
        } else {
          predbeta <- res0[[method]]$beta
        } 
        SE <- ((beta - predbeta)^2)
        AE <- abs(beta - predbeta)
        nz <- rowSums(exprs(e)) > 0
        de <- beta != 0
        RMSE <- sqrt(mean(SE[nz]))
        MAE <- mean(AE[nz])
        DiffRMSE <- sqrt(mean(SE[nz & de]))
        DiffMAE <- mean(AE[nz & de])
        data.frame(RMSE=RMSE, MAE=MAE, DiffRMSE=DiffRMSE, DiffMAE)
      }))
      
      data.frame(m=rep(m, length(methods)),
                 type=rep(type, length(methods)), 
                 method=methods, 
                 RMSE=resdf$RMSE,
                 MAE=resdf$MAE,
                 DiffRMSE=resdf$DiffRMSE,
                 DiffMAE=resdf$DiffMAE)
    }))
  }))
}))

res$method <- factor(res$method, methods)

save(res, file="results_simulateLFCAccuracy.RData")

