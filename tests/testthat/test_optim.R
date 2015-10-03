set.seed(1)
dds <- makeExampleDESeqDataSet(n=100)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
# make a large predictor to test scaling
colData(dds)$condition <- rnorm(ncol(dds),0,1000)
modelMatrix <- model.matrix(~ condition, as.data.frame(colData(dds)))
fit <- DESeq2:::fitNbinomGLMs(dds, modelMatrix=modelMatrix, 
                              modelFormula = ~ condition,
                              alpha_hat = dispersions(dds),
                              lambda = c(2,2),
                              renameCols=TRUE, betaTol=1e-8,
                              maxit=100, useOptim=TRUE,
                              useQR=TRUE, forceOptim=FALSE)
fitOptim <- DESeq2:::fitNbinomGLMs(dds, modelMatrix=modelMatrix, 
                                   modelFormula = ~ condition,
                                   alpha_hat = dispersions(dds),
                                   lambda = c(2,2),
                                   renameCols=TRUE, betaTol=1e-8,
                                   maxit=100, useOptim=TRUE,
                                   useQR=TRUE, forceOptim=TRUE)
expect_equal(fit$betaMatrix, fitOptim$betaMatrix,tolerance=1e-6)
expect_equal(fit$betaSE, fitOptim$betaSE,tolerance=1e-6)

# test optim gives same lfcSE
set.seed(1)
dds <- makeExampleDESeqDataSet(n=100, m=10)
counts(dds)[1,] <- c(rep(0L,5),c(1000L,1000L,0L,0L,0L))
dds <- DESeq(dds, betaPrior=FALSE)
res1 <- results(dds, contrast=c("condition","B","A"))
res2 <- results(dds, contrast=c(0,1))
expect_true(all.equal(res1$lfcSE, res2$lfcSE))
expect_true(all.equal(res1$pvalue, res2$pvalue))

