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

