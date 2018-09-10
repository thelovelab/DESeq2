context("tximport and tximeta")
test_that("tximport works", {
  library(tximport)
  library(readr)
  library(tximportData)
  dir <- system.file("extdata", package="tximportData")
  samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
  samples$condition <- factor(rep(c("A","B"),each=3))
  rownames(samples) <- samples$run
  samples[,c("pop","center","run","condition")]
  files <- file.path(dir,"salmon", samples$run, "quant.sf.gz")
  names(files) <- samples$run
  #tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
  txi <- tximport(files, type="salmon", txOut=TRUE)
  dds <- DESeqDataSetFromTximport(txi,
                                  colData = samples,
                                  design = ~ condition)
  expect_true("avgTxLength" %in% assayNames(dds))
  dds <- estimateSizeFactors(dds)
  expect_true("normalizationFactors" %in% assayNames(dds))
  txi2 <- txi
  # Note to users: this is NOT the ideal way to make CFA, instead use tximport()
  txi2$counts <- makeCountsFromAbundance(txi$counts,
                                         txi$abundance,
                                         txi$length,
                                         countsFromAbundance="lengthScaledTPM")
  txi2$countsFromAbundance <- "lengthScaledTPM"
  dds <- DESeqDataSetFromTximport(txi2,
                                  colData = samples,
                                  design = ~ condition)
  expect_true("counts" == assayNames(dds))
  dds <- estimateSizeFactors(dds)
  expect_true("counts" == assayNames(dds))
})
test_that("tximeta works", {
  library(tximportData)
  library(tximeta)
  dir <- system.file("extdata/salmon_dm", package="tximportData")
  files <- file.path(dir, "SRR1197474_cdna", "quant.sf.gz") 
  coldata <- data.frame(files, names="SRR1197474", condition="A", stringsAsFactors=FALSE)
  dir <- system.file("extdata", package="tximeta")
  indexDir <- file.path(dir, "Drosophila_melanogaster.BDGP6.cdna.v92_salmon_0.10.2")
  fastaFTP <- "ftp://ftp.ensembl.org/pub/release-92/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.cdna.all.fa.gz"
  dir2 <- system.file("extdata/salmon_dm", package="tximportData")
  gtfPath <- file.path(dir2,"Drosophila_melanogaster.BDGP6.92.gtf.gz")
  makeLinkedTxome(indexDir=indexDir,
                  source="Ensembl",
                  organism="Drosophila melanogaster",
                  release="92",
                  genome="BDGP6",
                  fasta=fastaFTP,
                  gtf=gtfPath,
                  write=FALSE)
  se <- tximeta(coldata)
  gse <- summarizeToGene(se)
  # warning about 1 file... ok
  expect_warning({ dds <- DESeqDataSet(gse, ~1) })
  expect_true("avgTxLength" %in% assayNames(dds))
  dds <- estimateSizeFactors(dds)
  expect_true("normalizationFactors" %in% assayNames(dds))
})
