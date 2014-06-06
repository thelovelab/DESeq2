test_results <- function() {
  # test some special cases for results()
  # using designs with a + 0 
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200,m=12)
  dds$condition <- factor(rep(1:3,each=4))
  dds$group <- factor(rep(1:2,length=ncol(dds)))

  counts(dds)[1,] <- rep(c(100L,200L,400L),each=4)
  
  design(dds) <- ~ condition + 0
  dds <- DESeq(dds, betaPrior=FALSE)

  checkEqualsNumeric(results(dds)[1,2], 2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","1"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","3","2"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","3"))[1,2], -2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","3"))[1,2], -1, tolerance=.1)
  
  design(dds) <- ~ group + condition + 0
  dds <- DESeq(dds, betaPrior=FALSE)

  checkEqualsNumeric(results(dds)[1,2], 2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","1"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","3","2"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","3"))[1,2], -2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","3"))[1,2], -1, tolerance=.1)
}
