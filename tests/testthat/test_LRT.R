context("LRT")
test_that("test='LRT' gives correct errors", {
  dds <- makeExampleDESeqDataSet(n=100, m=4)
  dds$group <- factor(c(1,2,1,2))
  design(dds) <- ~ condition
  expect_error(DESeq(dds, test="LRT", reduced=~group))
  expect_error(DESeq(dds, test="LRT", reduced=~1, modelMatrixType="expanded"))
  expect_error(DESeq(dds,test="LRT",reduced=~group, betaPrior=TRUE))
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  expect_error(nbinomLRT(dds))
})

test_that("glmGamPoi works", {

  dds <- makeExampleDESeqDataSet(n=100, m=20)
  dds$group <- factor(rep(1:2, times=10))
  design(dds) <- ~group + condition
  dds <- DESeq(dds, test="LRT", reduced=~1, fitType="glmGamPoi")

  # test outlier
  counts(dds)[1,1] <- 10000L
  design(dds) <- ~condition
  dds <- DESeq(dds, test="LRT", reduced=~1, fitType="glmGamPoi")
  
  # Michael Schubert's test
  dds <- makeExampleDESeqDataSet(n=100, m=8)
  nm = matrix(2, nrow=nrow(dds), ncol=ncol(dds))
  dds$group <- factor(rep(1:2,times=4))
  design(dds) <- ~group + condition
  dds <- estimateSizeFactors(dds, normMatrix=nm)
  dds <- estimateDispersions(dds, fitType="glmGamPoi")
  dds <- nbinomLRT(dds, reduced=~1, type="glmGamPoi")
  
})

test_that("test='LRT' with quasi-likelihood estimates gives correct errors", {
  
  dds <- makeExampleDESeqDataSet(n=100, m=4)
  dds$group <- factor(c(1,2,1,2))
  design(dds) <- ~ condition + group
  expect_warning(DESeq(dds, test = "Wald", fitType = "glmGamPoi"))
  dds <- estimateSizeFactors(dds)
  dds_gp <- estimateDispersions(dds)
  expect_error(nbinomLRT(dds_gp, reduced = ~ condition, type = "glmGamPoi"))
    
})
