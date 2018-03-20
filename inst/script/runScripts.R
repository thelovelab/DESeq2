runDESeq2 <- function(e, retDDS=FALSE) {
  counts <- exprs(e)
  mode(counts) <- "integer"
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds, quiet=TRUE)
  res <- results(dds)
  logFC <- res$log2FoldChange
  pval <- res$pvalue
  padj <- res$padj
  pval[is.na(pval)] <- 1
  pval[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pval=pval, padj=padj, logFC=logFC))
}

runDESeq2LRT <- function(e, retDDS=FALSE) {
  counts <- exprs(e)
  mode(counts) <- "integer"
  dds <- DESeqDataSetFromMatrix(counts, DataFrame(pData(e)), ~ condition)
  dds <- DESeq(dds,test="LRT",reduced=~1,quiet=TRUE)
  res <- results(dds)
  logFC <- res$log2FoldChange
  pval <- res$pvalue
  padj <- res$padj
  pval[is.na(pval)] <- 1
  pval[rowSums(exprs(e)) == 0] <- NA
  padj[is.na(padj)] <- 1
  return(list(pval=pval, padj=padj, logFC=logFC))
}

runEdgeR <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  ## dgel <- estimateGLMCommonDisp(dgel, design)
  ## dgel <- estimateGLMTrendedDisp(dgel, design)
  ## dgel <- estimateGLMTagwiseDisp(dgel, design)
  dgel <- estimateDisp(dgel, design)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predlogFC <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  predlogFC10 <- predFC(exprs(e), design, prior.count=10, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pval <- edger.lrt$table$PValue
  pval[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pval,method="BH")
  padj[is.na(padj)] <- 1
  list(pval=pval, padj=padj,
       logFC=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predlogFC=predlogFC[,"pData(e)$conditionB"],
       predlogFC10=predlogFC10[,"pData(e)$conditionB"])
}

runEdgeRRobust <- function(e) {
  design <- model.matrix(~ pData(e)$condition)
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  # settings for robust from robinson_lab/edgeR_robust/robust_simulation.R
  dgel <- estimateGLMRobustDisp(dgel, design, maxit=6)
  edger.fit <- glmFit(dgel, design)
  edger.lrt <- glmLRT(edger.fit)
  predlogFC <- predFC(exprs(e), design, offset=getOffset(dgel), dispersion=dgel$tagwise.dispersion)
  pval <- edger.lrt$table$PValue
  pval[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pval,method="BH")
  padj[is.na(padj)] <- 1
  list(pval=pval, padj=padj,
       logFC=log2(exp(1)) * edger.fit$coefficients[,"pData(e)$conditionB"],
       predlogFC=predlogFC[,"pData(e)$conditionB"])
}

runDSS <- function(e) {
  X <- as.matrix(exprs(e))
  colnames(X) <- NULL
  designs <- as.character(pData(e)$condition)
  seqData <- newSeqCountSet(X, designs)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  # typically gives locfdr warnings about df
  suppressWarnings({
    result <- waldTest(seqData, "B", "A")
  })
  result <- result[match(rownames(seqData),rownames(result)),]
  pval <- result$pval
  pval[rowSums(exprs(e)) == 0] <- NA
  # adjustment with BH (not the default FDR)
  padj <- p.adjust(pval,method="BH")
  padj[is.na(padj)] <- 1
  list(pval=pval, padj=padj, logFC=( log2(exp(1)) * result$lfc ))
}

runDSSFDR <- function(e) {
  X <- as.matrix(exprs(e))
  colnames(X) <- NULL
  designs <- as.character(pData(e)$condition)
  seqData <- newSeqCountSet(X, designs)
  seqData <- estNormFactors(seqData)
  seqData <- estDispersion(seqData)
  # typically gives locfdr warnings about df
  suppressWarnings({
    result <- waldTest(seqData, "B", "A")
  })
  result <- result[match(rownames(seqData),rownames(result)),]
  pval <- result$pval
  pval[rowSums(exprs(e)) == 0] <- NA
  padj <- result$fdr
  padj[is.na(padj)] <- 1
  list(pval=pval, padj=padj, logFC=( log2(exp(1)) * result$lfc ))
}

runVoom <- function(e) {
  design <- model.matrix(~ condition, pData(e))
  dgel <- DGEList(exprs(e))
  dgel <- calcNormFactors(dgel)
  v <- voom(dgel,design,plot=FALSE)
  fit <- lmFit(v,design)
  fit <- eBayes(fit)
  tt <- topTable(fit,coef=ncol(design),n=nrow(dgel),sort.by="none")
  pval <- tt$P.Value 
  pval[rowSums(exprs(e)) == 0] <- NA
  padj <- p.adjust(pval,method="BH")
  padj[is.na(padj)] <- 1
  list(pval=pval, padj=padj, logFC=tt$logFC)
}

runSAMseq <- function(e) {
  set.seed(1)
  x <- exprs(e)
  y <- pData(e)$condition
  capture.output({
    samfit <- SAMseq(x, y, resp.type = "Two class unpaired")
  })
  logFC <- log2(samfit$samr.obj$foldchange)
  pval <- samr.pvalues.from.perms(samfit$samr.obj$tt, samfit$samr.obj$ttstar)
  pval[rowSums(exprs(e)) == 0] <- NA
  # adjustment with BH (not the default FDR)
  padj <- p.adjust(pval,method="BH")
  padj[is.na(padj)] <- 1
  list(pval=pval,padj=padj,logFC=logFC)
}

runSAMseqFDR <- function(e) {
  set.seed(1)
  x <- exprs(e)
  y <- pData(e)$condition
  capture.output({
    samfit <- SAMseq(x, y, resp.type = "Two class unpaired", fdr.output=1)
  })
  padj <- rep(1,nrow(e))
  idx <- as.numeric(samfit$siggenes.table$genes.up[,"Gene Name"])
  padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.up[,"q-value(%)"])
  idx <- as.numeric(samfit$siggenes.table$genes.lo[,"Gene Name"])
  padj[idx] <- 1/100 * as.numeric(samfit$siggenes.table$genes.lo[,"q-value(%)"])
  logFC <- log2(samfit$samr.obj$foldchange)
  pval <- rep(1,nrow(e))
  list(pval=pval,padj=padj,logFC=logFC)
}

runEBSeq <- function(e) {
  sizes <- MedianNorm(exprs(e))
  out <- capture.output({
    suppressMessages({
      res <- EBTest(Data = exprs(e),
                    Conditions = pData(e)$condition,
                    sizeFactors = sizes,
                    maxround = 5)
    })
  })
  padj <- rep(1, nrow(exprs(e)))
  # we use 1 - PPDE for the FDR cutoff as this is recommended in the EBSeq vignette
  padj[match(rownames(res$PPMat), rownames(e))] <- res$PPMat[,"PPEE"]
  logFC <- rep(0, nrow(exprs(e)))
  pval <- rep(1,nrow(e))
  list(pval=pval, padj=padj, logFC=logFC)
}
