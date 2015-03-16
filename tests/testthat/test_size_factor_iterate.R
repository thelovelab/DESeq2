true.sf <- 2^(rep(c(-2,-1,0,0,1,2),each=2))
dds <- makeExampleDESeqDataSet(sizeFactors=true.sf, n=100)
cts <- counts(dds)
idx <- cbind(seq_len(nrow(cts)), sample(ncol(dds), nrow(cts), replace=TRUE))
cts[idx] <- 0L
cts[1,1] <- 1000000L # an outlier
counts(dds) <- cts
dds <- estimateSizeFactors(dds, type="iterate")
sf <- sizeFactors(dds)
coefs <- coef(lm(sf ~ true.sf))
expect_true(abs(coefs[1]) < .1)
expect_true(abs(coefs[2] - 1) < .1)

