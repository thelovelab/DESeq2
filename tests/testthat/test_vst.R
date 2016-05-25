dds <- makeExampleDESeqDataSet(n=100, m=4)
design(dds) <- ~ 1
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dds <- estimateDispersionsFit(dds, fitType="parametric")
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
dds <- estimateDispersionsFit(dds, fitType="local")
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
dds <- estimateDispersionsFit(dds, fitType="mean")
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)  

# test VST basics/errors
dds <- makeExampleDESeqDataSet(n=20, m=4)
colnames(dds) <- NULL
varianceStabilizingTransformation(dds)
head(varianceStabilizingTransformation(assay(dds)))
expect_error(getVarianceStabilizedData(dds))

# test just matrix
vsd <- varianceStabilizingTransformation(counts(dds))

# test fast VST based on subsampling
dds <- makeExampleDESeqDataSet(n=20000, m=10)
vsd <- vst(dds)
vsd <- vst(counts(dds))

# test VST and normalization factors
dds <- makeExampleDESeqDataSet(n=100, m=10, betaSD=1.5)
nf <- matrix(exp(rnorm(1000,0,.2)),ncol=10)
normalizationFactors(dds) <- nf
vsd <- varianceStabilizingTransformation(dds, fitType="local")
