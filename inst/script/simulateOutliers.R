source("makeSim.R")
load("meanDispPairs.RData")
library("Biobase")
library("DESeq2")
library("edgeR")
source("runScripts.R")
algos <- list("DESeq2"=runDESeq2Outliers,
              "edgeR"=runEdgeR,
              "edgeR-robust"=runEdgeRRobust)
namesAlgos <- names(algos)
methods <- c("DESeq2", "DESeq2-noFilt", "DESeq2-noRepl", "edgeR", "edgeR-robust")

set.seed(1)
padjVector <- seq(from=0, to=1, length=201)
pvalsVector <- seq(from=0, to=.4, length=201)
n <- 4000
percentOutliers <- seq(from=0,to=.15,length=4)
ms <- c(10,20)
nreps <- 10

res <- do.call(rbind, lapply(ms, function(m) {
  do.call(rbind, lapply(percentOutliers, function(pOut) {
    resList <- lapply(seq_len(nreps), function(i) {
      condition <- factor(rep(c("A","B"), each = m/2))
      x <- model.matrix(~ condition)
      beta <- c(rep(0, n * 8/10), sample(c(-1,1), n * 2/10, TRUE))
      mat <- makeSim(n,m,x,beta,meanDispPairs)$mat
      idx.i <- sample(n, round(n*pOut))
      idx.j <- sample(m, round(n*pOut), TRUE)
      mat[cbind(idx.i,idx.j)] <- mat[cbind(idx.i,idx.j)] * 100
      e <- ExpressionSet(mat, AnnotatedDataFrame(data.frame(condition)))
      resTest0 <- lapply(algos, function(f) f(e))
      # this avoids re-running DESeq2 without filtering or replacement
      resTest <- list()
      resTest[names(resTest0[["DESeq2"]])] <- resTest0[["DESeq2"]]
      resTest[c("edgeR","edgeR-robust")] <- resTest0[c("edgeR","edgeR-robust")]
      resTest[["beta"]] <- beta
      resTest[["nonzero"]] <- rowSums(exprs(e)) > 0
      resTest
    })
    sens <- sapply(methods, function(method) {
      sensMat <- sapply(seq_along(resList), function(i) {
        sapply(pvalsVector, function(p) {
          idx <- resList[[i]][["beta"]] != 0 & resList[[i]][["nonzero"]]
          mean((resList[[i]][[method]]$pvals < p)[idx])
        })
      })
      apply(sensMat, 1, median)
    })
    spec <- sapply(methods, function(method) {
      specMat <- sapply(seq_along(resList), function(i) {
        sapply(pvalsVector, function(p) {
          idx <- resList[[i]][["beta"]] == 0 & resList[[i]][["nonzero"]]
          mean((resList[[i]][[method]]$pvals >= p)[idx])
        })
      })
      apply(specMat, 1, median)
    })
    senspadj <- sapply(methods, function(method) {
      padjMat <- sapply(seq_along(resList), function(i) {
        sapply(pvalsVector, function(p) {
          idx <- resList[[i]][["nonzero"]]
          smallp <- resList[[i]][[method]]$pvals[idx] < p
          if (sum(smallp) == 0) 0 else max(resList[[i]][[method]]$padj[idx][ smallp ])
        })
      })
      apply(padjMat, 1, median)
    })   
    prec <- sapply(methods, function(method) {
      precMat <- sapply(seq_along(resList), function(i) {
        sapply(padjVector, function(p) {
          idx <- resList[[i]][[method]]$padj < p
          if (sum(idx) == 0) 1 else mean( (resList[[i]][["beta"]] != 0)[idx] )
        })
      })
      apply(precMat, 1, median)
    })
    data.frame(sensitivity = as.vector(sens),
               pvals = rep(pvalsVector, times=length(methods)),
               senspadj = as.vector(senspadj),
               oneminusspec = 1 - as.vector(spec),
               oneminusprec = 1 - as.vector(prec),
               precpadj = rep(padjVector, times=length(methods)),
               algorithm = rep(methods, each=length(pvalsVector)),
               m = rep(m, length(methods) * length(pvalsVector)),
               percentOutlier = rep(pOut, length(methods) * length(pvalsVector)))
  }))
}))

save(res, file="results_simulateOutliers.RData")

