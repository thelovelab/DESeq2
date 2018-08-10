context("construction_errors")
test_that("proper errors thrown in object construction", {
  coldata <- DataFrame(x=factor(c("A","A","B","B")),
                       xx=factor(c("A","A","B","B ")),
                       xwNA=factor(c("A","A","B",NA)),
                       name=letters[1:4],
                       ident=factor(rep("A",4)),
                       num=1:4,
                       missinglevels=factor(c("A","A","B","B"), levels=c("A","B","C")),
                       notref=factor(c("control","control","abc","abc")),
                       row.names=1:4)
  counts <- matrix(1:16, ncol=4)
  dup.rownms.counts <- matrix(1:16, ncol=4, dimnames=list(c(1,2,3,3),1:4))

  expect_message(DESeqDataSet(SummarizedExperiment(list(foo=counts), colData=coldata), ~ x), "renaming the first")
  expect_error(DESeqDataSetFromMatrix(matrix(c(1:11,-1),ncol=4), coldata, ~ x), "values in assay are negative")
  expect_error(DESeqDataSetFromMatrix(matrix(c(1:11,0.5),ncol=4), coldata, ~ x), "are not integers")
  expect_error(DESeqDataSetFromMatrix(matrix(rep(0,16),ncol=4), coldata, ~ x), "all samples have 0 counts")
  expect_warning(DESeqDataSetFromMatrix(matrix(rep(1:4,4),ncol=4), coldata, ~ x), "have equal values")
  expect_warning(DESeqDataSetFromMatrix(dup.rownms.counts, coldata, ~ x), "duplicate rownames")
  expect_error(DESeqDataSetFromMatrix(counts, colData=coldata, ~xwNA), "cannot contain NA")
  expect_error(DESeqDataSetFromMatrix(counts, coldata, ~ y), "must be columns in colData")
  expect_warning(DESeqDataSetFromMatrix(counts, coldata, ~ name), "are characters")
  expect_error(DESeqDataSetFromMatrix(counts, coldata, ~ ident), "all samples having the same value")
  expect_message(DESeqDataSetFromMatrix(counts, coldata, ~ num), "contains a numeric variable")
  expect_message(DESeqDataSetFromMatrix(counts, coldata, ~ missinglevels), "were dropped")
  expect_message(DESeqDataSetFromMatrix(counts, coldata, ~ notref), "not the reference level")
  expect_error(DESeqDataSetFromMatrix(counts, coldata, ~ident + x), "design contains")
  expect_message(DESeqDataSetFromMatrix(counts, coldata, ~xx), "characters other than")

  coldata2 <- data.frame(ord=ordered(rep(1:2,each=2)))
  expect_error(DESeqDataSetFromMatrix(counts, coldata2, ~ord), "ordered")
  
  chr.counts <- counts
  mode(chr.counts) <- "character"
  expect_error(DESeqDataSetFromMatrix(chr.counts, coldata, ~x), "should be numeric")
  
  # same colnames but in different order:
  expect_error(DESeqDataSetFromMatrix(matrix(1:16, ncol=4, dimnames=list(1:4, 4:1)), coldata, ~ x))

  # testing incoming metadata columns
  coldata <- DataFrame(x=factor(c("A","A","B","B")))
  rowranges <- GRanges("1", IRanges(1 + 0:3 * 10, width=10))
  se <- SummarizedExperiment(list(counts=counts), colData=coldata, rowRanges=rowranges)
  mcols(colData(se)) <- DataFrame(info="x is a factor")
  mcols(se)$id <- 1:4
  mcols(mcols(se)) <- DataFrame(info="the gene id")
  dds <- DESeqDataSet(se, ~ x)
  mcols(colData(dds))
  mcols(mcols(dds))

})
