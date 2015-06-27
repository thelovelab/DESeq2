# test disperion errors

dds <- makeExampleDESeqDataSet(n=100, m=2)
dds <- estimateSizeFactors(dds)
expect_error(estimateDispersionsGeneEst(dds))

set.seed(1)
dds <- makeExampleDESeqDataSet(n=100, m=4, dispMeanRel=function(x) 0.001 + x/1e3, interceptMean=8, interceptSD=2)
dds <- estimateSizeFactors(dds)
mcols(dds)$dispGeneEst <- rep(1e-7, 100)
expect_error(estimateDispersionsFit(dds))
dds <- estimateDispersionsGeneEst(dds)
expect_message(estimateDispersionsFit(dds))

dds <- makeExampleDESeqDataSet(n=100, m=4)
dds <- estimateSizeFactors(dds)
mcols(dds)$dispGeneEst <- rep(1e-7, 100)
dispersionFunction(dds) <- function(x) 1e-6
expect_warning(estimateDispersionsMAP(dds))

dds <- makeExampleDESeqDataSet(n=100, m=4)
dds <- estimateSizeFactors(dds)
levels(dds$condition) <- c("A","B","C")
expect_error(estimateDispersions(dds))
dds$condition <- droplevels(dds$condition)
dds$group <- dds$condition
design(dds) <- ~ group + condition
expect_error(estimateDispersions(dds))
