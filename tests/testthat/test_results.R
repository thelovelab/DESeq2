context("results")
test_that("results works as expected and throws errors", {
  ## test contrasts
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200,m=12)
  dds$condition <- factor(rep(1:3,each=4))
  dds$group <- factor(rep(1:2,length=ncol(dds)))
  counts(dds)[1,] <- rep(c(100L,200L,800L),each=4)

  design(dds) <- ~ group + condition

  # calling results too early
  expect_error(results(dds))

  sizeFactors(dds) <- rep(1, ncol(dds))
  dds <- DESeq(dds)
  head(coef(dds))
  res <- results(dds)
  show.res <- capture.output(show(res))
  summary.res <- capture.output(summary(res))

  # various results error checking
  expect_error(results(dds, test="LRT"))
  expect_error(results(dds, altHypothesis="lessAbs"))
  expect_error(results(dds, name=c("Intercept","group1")))
  expect_error(results(dds, contrast=c("foo","B","A")))
  expect_error(results(dds, contrast=c("condition","4","1")))
  expect_error(results(dds, test="foo"))
  expect_error(results(dds, contrast=FALSE))
  expect_error(results(dds, contrast=letters[1:4]))
  expect_error(results(dds, contrast=c("condition","1","1")))
  results(dds, independentFiltering=FALSE)
  results(dds, contrast=list("condition_2_vs_1"))
  expect_error(results(dds, contrast=list("condition_2_vs_1","condition_3_vs_1","condition_3_vs_1")))
  expect_error(results(dds, contrast=list("condition_2_vs_1",1)))
  expect_error(results(dds, contrast=list("condition_2_vs_1","foo")))
  expect_error(results(dds, contrast=list("condition_2_vs_1","condition_2_vs_1")))
  expect_error(results(dds, contrast=list(character(), character())))
  expect_error(results(dds, contrast=rep(0, 6)))

  # check to see if the contrasts with expanded model matrix
  # are close to expected (although shrunk due to the beta prior).
  # lfcShrink() here calls results()
  lfc31 <- lfcShrink(dds,contrast=c("condition","3","1"))[1,"log2FoldChange"]
  lfc21 <- lfcShrink(dds,contrast=c("condition","2","1"))[1,"log2FoldChange"]
  lfc32 <- lfcShrink(dds,contrast=c("condition","3","2"))[1,"log2FoldChange"]
  expect_equal(lfc31, 3, tolerance=.1)
  expect_equal(lfc21, 1, tolerance=.1)
  expect_equal(lfc32, 2, tolerance=.1)
  expect_equal(results(dds,contrast=c("condition","1","3"))[1,2], -3, tolerance=.1)
  expect_equal(results(dds,contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  expect_equal(results(dds,contrast=c("condition","2","3"))[1,2], -2, tolerance=.1)

  # check that results are not changed by releveling
  dds2 <- dds
  colData(dds2)$condition <- relevel(colData(dds2)$condition, "2")
  dds2 <- DESeq(dds2)
  expect_equal(lfcShrink(dds2,contrast=c("condition","3","1"))[1,"log2FoldChange"],
               lfc31, tolerance=1e-6)
  expect_equal(lfcShrink(dds2,contrast=c("condition","2","1"))[1,"log2FoldChange"],
               lfc21, tolerance=1e-6)
  expect_equal(lfcShrink(dds2,contrast=c("condition","3","2"))[1,"log2FoldChange"],
               lfc32, tolerance=1e-6)

  # test a number of contrast as list options
  expect_equal(results(dds, contrast=list("condition_3_vs_1","condition_2_vs_1"))[1,2],
               2, tolerance=1e-6)
  results(dds, contrast=list("condition_3_vs_1","condition_2_vs_1"), listValues=c(.5,-.5))
  results(dds, contrast=list("condition_3_vs_1",character()))
  results(dds, contrast=list("condition_3_vs_1",character()), listValues=c(.5,-.5))
  results(dds, contrast=list(character(),"condition_2_vs_1"))
  results(dds, contrast=list(character(),"condition_2_vs_1"), listValues=c(.5,-.5))

  # test no prior on intercept
  expect_equivalent(attr(dds,"betaPriorVar"), rep(1e6, 4))

  # test thresholding
  resLFC <- results(dds, lfcThreshold=log2(1.5))
  results(dds, lfcThreshold=1, altHypothesis="lessAbs")
  results(dds, lfcThreshold=1, altHypothesis="greater")
  results(dds, lfcThreshold=1, altHypothesis="less")

  summary.res <- capture.output({ summary(resLFC)})
  expect_true(any(grepl("0.58", summary.res)))
  
  dds3 <- DESeq(dds, betaPrior=TRUE)
  expect_error(results(dds3, lfcThreshold=1, altHypothesis="lessAbs"))
})

test_that("results: designs with zero intercept", {
  # test some special cases for results()
  # using designs with +0 
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100,m=12)
  dds$condition <- factor(rep(1:3,each=4))
  dds$group <- factor(rep(1:2,length=ncol(dds)))

  counts(dds)[1,] <- rep(c(100L,200L,400L),each=4)

  design(dds) <- ~ condition + 0
  dds <- DESeq(dds, betaPrior=FALSE)

  expect_equal(results(dds)[1,2], 2, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","2","1"))[1,2], 1, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","3","2"))[1,2], 1, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","1","3"))[1,2], -2, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","2","3"))[1,2], -1, tolerance=.1)
  expect_error(results(dds, contrast=c("condition","4","1")))

  design(dds) <- ~ group + condition + 0
  dds <- DESeq(dds, betaPrior=FALSE)

  expect_equal(results(dds)[1,2], 2, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","2","1"))[1,2], 1, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","3","2"))[1,2], 1, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","1","3"))[1,2], -2, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  expect_equal(results(dds, contrast=c("condition","2","3"))[1,2], -1, tolerance=.1)
})

test_that("results: likelihood ratio test", {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100)
  dds$group <- factor(rep(1:2,6))
  design(dds) <- ~ group + condition
  dds <- DESeq(dds, test="LRT", reduced=~group)

  expect_true(!all(results(dds,name="condition_B_vs_A")$stat ==
              results(dds,name="condition_B_vs_A",test="Wald")$stat))

  # LFC are already MLE
  expect_error(results(dds, addMLE=TRUE))
  expect_error(results(dds, lfcThreshold=1, test="LRT"))

  expect_true(all(results(dds, test="LRT", contrast=c("group","1","2"))$log2FoldChange ==
              -1 * results(dds, test="LRT", contrast=c("group","2","1"))$log2FoldChange))
})

test_that("results basics regarding format, tidy, MLE, remove are working", {
  dds <- makeExampleDESeqDataSet(n=100)
  dds <- DESeq(dds)
  res <- results(dds, format="GRanges")
  expect_warning(results(dds, format="GRangesList"))

  rowRanges(dds) <- as(rowRanges(dds), "GRangesList")
  dds <- DESeq(dds)
  expect_message(results(dds, format="GRanges"))

  # check tidy-ness
  res <- results(dds, tidy=TRUE)
  expect_true(colnames(res)[1] == "row")
  expect_true(is(res, "data.frame"))

  # test MLE and 'name'
  dds2 <- DESeq(dds, betaPrior=TRUE)
  results(dds2, addMLE=TRUE)
  expect_error(results(dds, name="condition_B_vs_A", addMLE=TRUE))

  # test remove results
  dds <- removeResults(dds)
  expect_true(!any(mcols(mcols(dds))$type == "results"))
})

test_that("custom filters can be provided to results()", {
  # try a custom filter function
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200, m=4, betaSD=rep(c(0,2),c(150,50)))
  dds <- DESeq(dds)
  res <- results(dds)
  method <- "BH"
  alpha <- 0.1

  customFilt <- function(res, filter, alpha, method) {
    if (missing(filter)) {
      filter <- res$baseMean
    }
    theta <- 0:10/10
    cutoff <- quantile(filter, theta)
    numRej <- sapply(cutoff, function(x) sum(p.adjust(res$pvalue[filter > x]) < alpha, na.rm=TRUE))
    threshold <- theta[which(numRej == max(numRej))[1]]
    res$padj <- numeric(nrow(res))
    idx <- filter > quantile(filter, threshold)
    res$padj[!idx] <- NA
    res$padj[idx] <- p.adjust(res$pvalue[idx], method=method)
    res
  }

  resCustom <- results(dds, filterFun=customFilt)
  #plot(res$padj, resCustom$padj);abline(0,1)
})
