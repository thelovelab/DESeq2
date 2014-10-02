source("makeSim.R")
load("meanDispPairs.RData")
library("DESeq2")
library("PoiClaClu")
library("mclust")

set.seed(1)
n <- 2000
# create 20 samples, then remove first group, leaving 16
# this way the groups are equidistant from each other
m <- 20
k <- 4
methods <- c("norm Eucl","log2 Eucl","rlog Eucl","VST Eucl","PoisDist")
condition0 <- factor(rep(c("null","A","B","C","D"), each = m/(k+1)))
x <- model.matrix(~ condition0)
rnormsds <- list(seq(from=0, to=.6, length=7),
                 seq(from=0, to=.8, length=7),
                 seq(from=0, to=1.2, length=7))
sfs <- list(equal=rep(1,m), unequal=rep(c(1,1,1/3,3), times=(k+1)))
dispScales <- c(.1, .25, 1)
nreps <- 20

res <- do.call(rbind, lapply(seq_along(dispScales), function(idx) {
  dispScale <- dispScales[idx]
  do.call(rbind, lapply(rnormsds[[idx]], function(rnormsd) {
    do.call(rbind, lapply(seq_along(sfs), function(sf.idx) {
      sf <- sfs[[sf.idx]]
      do.call(rbind, lapply(seq_len(nreps), function(i) {
        beta <- replicate(k, c(rep(0,8/10 * n), rnorm(2/10 * n, 0, rnormsd)))
        mdp <- meanDispPairs
        mdp$disp <- mdp$disp * dispScale
        mat0 <- makeSim(n,m,x,beta,mdp,sf)$mat
        mat <- mat0[,5:20]
        mode(mat) <- "integer"
        condition <- droplevels(condition0[5:20])
        dds <- DESeqDataSetFromMatrix(mat, DataFrame(condition), ~ 1)
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersionsGeneEst(dds)
        # don't warn if local fit is used
        dds <- suppressWarnings({estimateDispersionsFit(dds)})
        norm <- t(counts(dds, normalized=TRUE))
        lognorm <- t(log2(counts(dds, normalized=TRUE) + 1))
        rld <- t(assay(rlog(dds, blind=FALSE)))
        vsd <- t(assay(varianceStabilizingTransformation(dds, blind=FALSE)))
        poiDist <- PoissonDistance(t(mat))$dd

        normARI <- adjustedRandIndex(condition, cutree(hclust(dist(norm)),k=k))
        lognormARI <- adjustedRandIndex(condition, cutree(hclust(dist(lognorm)),k=k))
        rlogARI <- adjustedRandIndex(condition, cutree(hclust(dist(rld)),k=k))
        vstARI <- adjustedRandIndex(condition, cutree(hclust(dist(vsd)),k=k))
        poiDistARI <- adjustedRandIndex(condition, cutree(hclust(poiDist),k=k))

        data.frame(ARI = c(normARI, lognormARI, rlogARI, vstARI, poiDistARI),
                   method = methods,
                   rnormsd = rep(rnormsd,length(methods)),
                   dispScale = rep(dispScale,length(methods)),
                   sizeFactor = rep(names(sfs)[sf.idx], length(methods)))
      }))
    }))
  }))
}))

res$method <- factor(res$method, methods)

save(res, file="results_simulateCluster.RData")

