source("makeSim.R")
load("meanDispPairs.RData")
library("Biobase")
library("DESeq2")
library("edgeR")
library("limma")
library("samr")
library("DSS")
library("EBSeq")
source("runScripts.R")
algos <- list("DESeq2"=runDESeq2,
              "edgeR"=runEdgeR,"edgeR-robust"=runEdgeRRobust,
              "DSS"=runDSS,"DSS-FDR"=runDSSFDR,
              "voom"=runVoom,
              "SAMseq"=runSAMseq,"SAMseq-FDR"=runSAMseqFDR,
              "EBSeq"=runEBSeq)
namesAlgos <- names(algos)

n <- 10000
effSizeLevels <- log2(c(2,3,4))
mLevels <- c(6,8,10,20)
nreps <- 6

effSizes <- rep(rep(effSizeLevels, each=nreps), times=length(mLevels))
ms <- rep(mLevels, each=nreps * length(effSizeLevels))

library("BiocParallel")
register(SerialParam())
#register(MulticoreParam(workers=8,verbose=TRUE))

resList <- bplapply(seq_along(ms), function(i) {
  set.seed(i)
  m <- ms[i]
  es <- effSizes[i]
  condition <- factor(rep(c("A","B"), each = m/2))
  x <- model.matrix(~ condition)
  beta <- c(rep(0, n * 8/10), sample(c(-es,es), n * 2/10, TRUE))
  mat <- makeSim(n,m,x,beta,meanDispPairs)$mat
  e <- ExpressionSet(mat, AnnotatedDataFrame(data.frame(condition)))
  resTest <- lapply(algos, function(f) f(e))
  nonzero <- rowSums(exprs(e)) > 0
  sensidx <- abs(beta) > 0 & nonzero
  sens <- sapply(resTest, function(z) mean((z$padj < .1)[sensidx]))
  rmf <- cut(rowMeans(mat), c(0, 20, 100, 300, Inf), include.lowest=TRUE)
  levels(rmf) <- paste0("sens",c("0to20","20to100","100to300","more300"))
  sensStratified <- t(sapply(resTest, function(z) tapply((z$padj < .1)[sensidx], rmf[sensidx], mean)))
  oneminusspecpvals <- sapply(resTest, function(z) mean((z$pvals < .01)[beta == 0 & nonzero], na.rm=TRUE))
  oneminusspecpadj <- sapply(resTest, function(z) mean((z$padj < .1)[beta == 0 & nonzero], na.rm=TRUE))
  oneminusprec <- sapply(resTest, function(z) {
      idx <- which(z$padj < .1)
      ifelse(sum(idx) == 0, 0, mean((beta == 0)[idx]))
  })
  data.frame(sensitivity=sens,
             sensStratified,
             oneminusspecpvals=oneminusspecpvals,
             oneminusspecpadj=oneminusspecpadj,
             oneminusprec=oneminusprec,
             algorithm=namesAlgos, effSize=es, m=m)
})
res <- do.call(rbind, resList)

res$algorithm <- factor(res$algorithm, namesAlgos)

sessInfo <- sessionInfo()

save(res, namesAlgos, sessInfo, file="results_simulateDE.RData")

