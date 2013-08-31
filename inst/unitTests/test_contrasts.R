test_contrasts <- function() {
  set.seed(1)
  n <- 1000
  dds <- makeExampleDESeqDataSet(n=n)
  counts(dds) <- matrix(as.integer(rpois(n*12,
                        lambda=rep(c(100,200,800),each=n*4))),ncol=12)
  colData(dds)$condition <- factor(rep(c("A","B","C"),each=4))
  sizeFactors(dds) <- rep(1,12)
  dispersions(dds) <- .1
  dds <- nbinomWaldTest(dds,betaPrior=FALSE)
  lfcCA <- results(dds)$log2FoldChange
  lfcBA <- results(dds,name="condition_B_vs_A")$log2FoldChange
  lfcCB <- results(dds,contrast=c("condition","C","B"))$log2FoldChange
  checkEqualsNumeric(median(lfcCA), 3, tolerance=.01)
  checkEqualsNumeric(median(lfcBA), 1, tolerance=.01)
  checkEqualsNumeric(median(lfcCB), 2, tolerance=.01)
  
  dds <- nbinomWaldTest(dds,betaPrior=FALSE,useQR=FALSE)
  lfcCA <- results(dds)$log2FoldChange
  lfcBA <- results(dds,name="condition_B_vs_A")$log2FoldChange
  lfcCB <- results(dds,contrast=c("condition","C","B"))$log2FoldChange
  checkEqualsNumeric(median(lfcCA), 3, tolerance=.01)
  checkEqualsNumeric(median(lfcBA), 1, tolerance=.01)
  checkEqualsNumeric(median(lfcCB), 2, tolerance=.01)
}
