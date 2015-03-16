set.seed(1)
dds <- makeExampleDESeqDataSet(n=100)
design(dds) <- ~ 1
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

ddsNew <- makeExampleDESeqDataSet(m=1,n=100)
counts(ddsNew)[,1] <- counts(dds)[,1]
sizeFactors(ddsNew)[1] <- sizeFactors(dds)[1]

# VST
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
dispersionFunction(ddsNew) <- dispersionFunction(dds)
vsdNew <- varianceStabilizingTransformation(ddsNew, blind=FALSE)
expect_equal(assay(vsd)[,1],assay(vsdNew)[,1],tolerance=1e-3)

# rlog
rld <- rlogTransformation(dds, blind=FALSE)  
mcols(ddsNew)$dispFit <- mcols(dds)$dispFit
betaPriorVar <- attr(rld,"betaPriorVar")
intercept <- mcols(rld)$rlogIntercept
rldNew <- rlogTransformation(ddsNew, blind=FALSE,
                             betaPriorVar=betaPriorVar,
                             intercept=intercept)
expect_equal(assay(rld)[,1],assay(rldNew)[,1],tolerance=1e-3)

# rlog fast
rld <- rlogTransformation(dds, blind=FALSE, fast=TRUE)  
mcols(ddsNew)$dispFit <- mcols(dds)$dispFit
B <- attr(rld,"B")
intercept <- mcols(rld)$rlogIntercept
rldNew <- rlogTransformation(ddsNew, blind=FALSE,
                             intercept=intercept, B=B,
                             fast=TRUE)
expect_equal(assay(rld)[,1],assay(rldNew)[,1],tolerance=1e-3)

