test_contrasts <- function() {
  set.seed(1)
  n <- 1000
  dds <- makeExampleDESeqDataSet(n=n)
  counts(dds) <- matrix(as.integer(rpois(n*12,
                        lambda=rep(c(100,200,800),each=n*4))),ncol=12)
  colData(dds)$condition <- factor(rep(c("A","B","C"),each=4))
  colData(dds)$group <- factor(rep(c("x","x","y","y"),3))
  design(dds) <- ~ group + condition
  sizeFactors(dds) <- rep(1,12)
  dispersions(dds) <- .1

  # check to see if the contrasts with expanded model matrix
  # are close to expected (although shrunk due to the beta prior)
  dds <- nbinomWaldTest(dds)
  lfcCA <- results(dds,contrast=c("condition","C","A"))$log2FoldChange
  lfcBA <- results(dds,contrast=c("condition","B","A"))$log2FoldChange
  lfcCB <- results(dds,contrast=c("condition","C","B"))$log2FoldChange
  checkEqualsNumeric(median(lfcCA), 3, tolerance=.1)
  checkEqualsNumeric(median(lfcBA), 1, tolerance=.1)
  checkEqualsNumeric(median(lfcCB), 2, tolerance=.1)

  # check that results are not changed by releveling
  dds2 <- dds
  colData(dds2)$condition <- relevel(colData(dds2)$condition, "B")
  dds2 <- nbinomWaldTest(dds2)
  lfcCA2 <- results(dds2,contrast=c("condition","C","A"))$log2FoldChange
  lfcBA2 <- results(dds2,contrast=c("condition","B","A"))$log2FoldChange
  lfcCB2 <- results(dds2,contrast=c("condition","C","B"))$log2FoldChange  
  checkEqualsNumeric(median(lfcCA), median(lfcCA2), tolerance=1e-6)
  checkEqualsNumeric(median(lfcBA), median(lfcBA2), tolerance=1e-6)
  checkEqualsNumeric(median(lfcCB), median(lfcCB2), tolerance=1e-6)

  design(dds) <- ~ group + condition + condition:group
  dds <- nbinomWaldTest(dds)
  # check the default prior  variance on the intercept and group LFC's
  checkEquals(attr(dds,"betaPriorVar")[1:6],
              c(Intercept=1e6, groupx=1e3, groupy=1e3,
                conditionA=1e3, conditionB=1e3, conditionC=1e3))
}
