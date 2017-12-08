context("tximport")
test_that("tximport works", {
  library("tximport")
  library("tximportData")
  library("readr")
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
  names(files) <- paste0("sample",1:6)
  tx2gene <- read.csv(file.path(dir, "tx2gene.gencode.v27.csv"))
  
  txi <- tximport(files, type="salmon", tx2gene=tx2gene)

  dds <- DESeqDataSetFromTximport(txi, samples, ~1)

  # test library size correction taking into account
  # the average transcript lengths
  dds <- estimateSizeFactors(dds)
  
  # test fpkm
  exprs <- fpm(dds)
  exprs <- fpkm(dds)

  # test length of 0
  txi2 <- txi
  txi2$length[1,1] <- 0
  expect_error(dds2 <- DESeqDataSetFromTximport(txi2, samples, ~1), "lengths")
})
