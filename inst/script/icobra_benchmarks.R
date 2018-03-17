suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(samr))
suppressPackageStartupMessages(library(DSS))
suppressPackageStartupMessages(library(EBSeq))
suppressPackageStartupMessages(library(iCOBRA))

source("makeSim.R")
source("runScripts.R")
algos <- list("DESeq2"=runDESeq2,
              "edgeR"=runEdgeR,
              "DSS"=runDSSFDR,
              "voom"=runVoom,
              "SAMseq"=runSAMseqFDR,
              "EBSeq"=runEBSeq)

n <- 5000
mLevels <- c(6, 8, 10) # total sample size
derLevels <- c(.1, .2, .3) # DE ratio

grid <- expand.grid(
  m=mLevels,
  der=derLevels)

library("parallel")
options(mc.cores=4)

# only simulate with base mean > 10
load("meanDispPairs_bottomly.rda")
mdp <- meanDispPairs[meanDispPairs$mean > 10,]

cdList <- mclapply(seq_len(nrow(grid)), function(i) {
  set.seed(i)
  m <- grid$m[i]
  der <- grid$der[i]
  condition <- factor(rep(c("A","B"), each = m/2))
  x <- model.matrix(~ condition)
  beta <- c(rep(0, (1 - der) * n),
            sample(c(-1,1),der * n,TRUE) * log2(runif(der * n,1.5,2)))
  status <- as.integer(beta != 0)
  truth <- data.frame(status=status, logFC=beta)
  mat <- makeSim(n,m,x,beta,mdp)$mat
  e <- ExpressionSet(mat, AnnotatedDataFrame(data.frame(condition)))
  system.time({
    resTest <- lapply(algos, function(f) f(e))
  })
  cd <- COBRAData(pval = data.frame(lapply(resTest, `[[`, "pval")),
                  padj = data.frame(lapply(resTest, `[[`, "padj")),
                  score = data.frame(lapply(resTest, `[[`, "logFC")),
                  truth = truth)
  cd@score$EBSeq <- cd@score$edgeR # inject LFC estimate for EBSeq
  cd
})

plotList <- list()
for (i in seq_len(nrow(grid))) {
  cat(i, "")
  cp <- calculate_performance(cdList[[i]],
                              binary_truth = "status", cont_truth = "logFC",
                              thrs=c(.05, .1),
                              aspects=c("fdrtpr","fdrtprcurve"))
  cobraplot <- prepare_data_for_plot(cp)
  title <- paste0("m=", grid$m[i], ";  ratio=", grid$der[i])
  plotList[[i]] <- plot_fdrtprcurve(cobraplot, title=title,
                                    xaxisrange=c(0, 0.3),
                                    yaxisrange=c(0.5, 1))
}

library(cowplot)
png(file="icobra.png", width=1200, height=1200)
do.call(plot_grid, plotList)
dev.off()

sessInfo <- sessionInfo()
save(cdList, sessInfo, file="icobra.rda")

if (FALSE) {
  load("bottomly_sumexp.RData")
  bottomly.se <- updateObject(bottomly)
  levels(bottomly.se$strain) <- c("C","D")
  bottomly.se$batch <- factor(bottomly.se$experiment.number)
  bottomly.se <- DESeqDataSet(bottomly.se, ~batch + strain)
  bottomly.se <- estimateSizeFactors(bottomly.se)
  bottomly.se <- estimateDispersions(bottomly.se)
  meanDispPairs <- with(mcols(bottomly.se),
                        data.frame(mean=baseMean, disp=dispGeneEst))
  meanDispPairs <- meanDispPairs[meanDispPairs$disp > 1e-6,]
  meanDispPairs <- meanDispPairs[!is.na(meanDispPairs$mean),]
  save(meanDispPairs, file="meanDispPairs_bottomly.rda")
}
