#' Apply a 'regularized log' transformation
#'
#' This function transforms the count data to the log2 scale in a way 
#' which minimizes differences between samples for rows with small counts,
#' and which normalizes with respect to library size.
#' The rlog transformation produces a similar variance stabilizing effect as
#' \code{\link{varianceStabilizingTransformation}},
#' though \code{rlog} is more robust in the
#' case when the size factors vary widely.
#' The transformation is useful when checking for outliers
#' or as input for machine learning techniques
#' such as clustering or linear discriminant analysis.
#' \code{rlog} takes as input a \code{\link{DESeqDataSet}} and returns a
#' \code{\link{SummarizedExperiment}} object.
#'
#' Note that neither rlog transformation nor the VST are used by the
#' differential expression estimation in \code{\link{DESeq}}, which always
#' occurs on the raw count data, through generalized linear modeling which
#' incorporates knowledge of the variance-mean dependence. The rlog transformation
#' and VST are offered as separate functionality which can be used for visualization,
#' clustering or other machine learning tasks. See the transformation section of the
#' vignette for more details.
#'
#' The transformation does not require that one has already estimated size factors
#' and dispersions.
#'
#' The regularization is on the log fold changes of the count for each sample
#' over an intercept, for each gene. As nearby count values for low counts genes
#' are almost as likely as the observed count, the rlog shrinkage is greater for low counts.
#' For high counts, the rlog shrinkage has a much weaker effect.
#' The fitted dispersions are used rather than the MAP dispersions
#' (so similar to the \code{\link{varianceStabilizingTransformation}}).
#' 
#' The prior variance for the shrinkag of log fold changes is calculated as follows: 
#' a matrix is constructed of the logarithm of the counts plus a pseudocount of 0.5,
#' the log of the row means is then subtracted, leaving an estimate of
#' the log fold changes per sample over the fitted value using only an intercept.
#' The prior variance is then calculated by matching the upper quantiles of the observed 
#' log fold change estimates with an upper quantile of the normal distribution.
#' A GLM fit is then calculated using this prior. It is also possible to supply the variance of the prior.
#' See the vignette for an example of the use and a comparison with \code{varianceStabilizingTransformation}.
#'
#' The transformed values, rlog(K), are equal to
#' \eqn{rlog(K_{ij}) = \log_2(q_{ij}) = \beta_{i0} + \beta_{ij}}{rlog(K_ij) = log2(q_ij) = beta_i0 + beta_ij},
#' with formula terms defined in \code{\link{DESeq}}.
#'
#' The parameters of the rlog transformation from a previous dataset
#' can be frozen and reapplied to new samples. See the 'Data quality assessment'
#' section of the vignette for strategies to see if new samples are
#' sufficiently similar to previous datasets. 
#' The frozen rlog is accomplished by saving the dispersion function,
#' beta prior variance and the intercept from a previous dataset,
#' and running \code{rlog} with 'blind' set to FALSE
#' (see example below).
#' 
#' @aliases rlog rlogTransformation
#' @rdname rlog
#' @name rlog
#' 
#' @param object a DESeqDataSet, or matrix of counts
#' @param blind logical, whether to blind the transformation to the experimental
#' design. blind=TRUE should be used for comparing samples in an manner unbiased by
#' prior information on samples, for example to perform sample QA (quality assurance).
#' blind=FALSE should be used for transforming data for downstream analysis,
#' where the full use of the design information should be made.
#' blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated.
#' If many of genes have large differences in counts due to
#' the experimental design, it is important to set blind=FALSE for downstream
#' analysis. 
#' @param fast if TRUE, an alternative rlog-like transformation is made,
#' wherein an optimal amount of linear shrinkage of sample-wise log fold changes
#' is calculated for each gene. 'Optimal' is defined by maximizing the same
#' posterior as in the standard rlog transformation.
#' @param intercept by default, this is not provided and calculated automatically.
#' if provided, this should be a vector as long as the number of rows of object,
#' which is log2 of the mean normalized counts from a previous dataset.
#' this will enforce the intercept for the GLM, allowing for a "frozen" rlog
#' transformation based on a previous dataset.
#' @param betaPriorVar a single value, the variance of the prior on the sample
#' betas, which if missing is estimated from the data
#' @param B for the \code{fast} method only.
#' by default, this is not provided and calculated automatically.
#' if provided, this should be a vector as long as the number of rows of object.
#' this is the amount of shrinkage from 0 to 1 for each row, and takes precedence
#' over internal calculation of B using betaPriorVar.
#' @param fitType in case dispersions have not yet been estimated for \code{object},
#' this parameter is passed on to \code{\link{estimateDispersions}} (options described there).
#' 
#' @return a \code{\link{DESeqTransform}} if a \code{DESeqDataSet} was provided,
#' or a matrix if a count matrix was provided as input.
#' Note that for \code{\link{DESeqTransform}} output, the matrix of
#' transformed values is stored in \code{assay(rld)}.
#' To avoid returning matrices with NA values, in the case of a row
#' of all zeros, the rlog transformation returns zeros
#' (essentially adding a pseudocount of 1 only to these rows).
#'
#' @references
#'
#' Reference for regularized logarithm (rlog):
#' 
#' Michael I Love, Wolfgang Huber, Simon Anders: Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 2014, 15:550. \url{http://dx.doi.org/10.1186/s13059-014-0550-8}
#' 
#' @seealso \code{\link{plotPCA}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=6,betaSD=1)
#' rld <- rlog(dds)
#' dists <- dist(t(assay(rld)))
#' plot(hclust(dists))
#'
#' # run the rlog transformation on one dataset
#' design(dds) <- ~ 1
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' rld <- rlog(dds, blind=FALSE)
#'
#' # apply the parameters to a new sample
#' 
#' ddsNew <- makeExampleDESeqDataSet(m=1)
#' mcols(ddsNew)$dispFit <- mcols(dds)$dispFit
#' betaPriorVar <- attr(rld,"betaPriorVar")
#' intercept <- mcols(rld)$rlogIntercept
#' rldNew <- rlog(ddsNew, blind=FALSE,
#'                intercept=intercept,
#'                betaPriorVar=betaPriorVar)
#'                            
#' 
#' @export
rlog <- function(object, blind=TRUE, fast=FALSE,
                 intercept, betaPriorVar, B, fitType="parametric") {
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names=colnames(object)), ~ 1)
  } else {
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~ 1
  }
  # sparsity test
  if (missing(intercept) & missing(B)) {
    sparseTest(counts(object, normalized=TRUE), .9, 100, .1)
  }
  if (blind | is.null(mcols(object)$dispFit)) {
    # estimate the dispersions on all genes, or if fast=TRUE subset to 1000 non-zero genes
    if (is.null(mcols(object)$baseMean)) {
      object <- getBaseMeansAndVariances(object)
    }
    if (!fast | sum(mcols(object)$baseMean > 0) <= 1000) {
      object <- estimateDispersionsGeneEst(object, quiet=TRUE)
      object <- estimateDispersionsFit(object, fitType, quiet=TRUE)
    } else {
      # select 1000 genes along the range of non-zero base means
      idx <- order(mcols(object)$baseMean)[round(seq(from=(sum(mcols(object)$baseMean==0) + 1),
                                                     to=nrow(object), length=1000))]      
      objectSub <- object[idx,]
      objectSub <- estimateDispersionsGeneEst(objectSub, quiet=TRUE)
      objectSub <- estimateDispersionsFit(objectSub, fitType, quiet=TRUE)
      # fill in the fitted dispersions for all genes
      mcols(object)$dispFit <- dispersionFunction(objectSub)(mcols(object)$baseMean)
    }
  }
  if (!missing(intercept)) {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
  }
  if (!missing(B)) {
    if (length(B) != nrow(object)) {
      stop("B should be as long as the number of rows of object")
    }
    if (!all(B >= 0 & B <= 1)) {
      stop("B should be defined between 0 and 1")
    }
  }
  if (fast) {
    rld <- rlogDataFast(object, intercept, betaPriorVar, B)
  } else {
    rld <- rlogData(object, intercept, betaPriorVar)
  }
  if (matrixIn) {
    return(rld)
  }
  se <- SummarizedExperiment(
           assays = rld,
           colData = colData(object),
           rowRanges = rowRanges(object),
           exptData = exptData(object))
  dt <- DESeqTransform(se)
  attr(dt,"betaPriorVar") <- attr(rld, "betaPriorVar")
  if (!is.null(attr(rld,"intercept"))) {
    mcols(dt)$rlogIntercept <- attr(rld,"intercept")
  }
  if (!is.null(attr(rld,"B"))) {
    attr(dt,"B") <- attr(rld,"B")
  }
  dt
}

#' @rdname rlog
#' @export
rlogTransformation <- rlog

###################### unexported


rlogData <- function(object, intercept, betaPriorVar) {
  if (is.null(mcols(object)$dispFit)) {
    stop("first estimate dispersion")
  }
  samplesVector <- as.character(seq_len(ncol(object)))
  if (!missing(intercept)) {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
  }
  if (is.null(mcols(object)$allZero) | is.null(mcols(object)$baseMean)) {
    object <- getBaseMeansAndVariances(object)
  }
  
  # make a design matrix with a term for every sample
  # this would typically produce unidentifiable solution
  # for the GLM, but we add priors for all terms except
  # the intercept
  samplesVector <- factor(samplesVector,levels=unique(samplesVector))
  if (missing(intercept)) {
    samples <- factor(c("null_level",as.character(samplesVector)),
                      levels=c("null_level",levels(samplesVector)))
    modelMatrix <- model.matrix(~samples)[-1,]
    modelMatrixNames <- colnames(modelMatrix)
    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  } else {
    # or we want to set the intercept using the
    # provided intercept instead
    samples <- factor(samplesVector)
    if (length(samples) > 1) {
      modelMatrix <- model.matrix(~ 0 + samples)
    } else {
      modelMatrix <- matrix(1,ncol=1)
      modelMatrixNames <- "samples1"
    }
    modelMatrixNames <- colnames(modelMatrix)
    if (!is.null(normalizationFactors(object))) { 
      nf <- normalizationFactors(object)
    } else {
      sf <- sizeFactors(object)
      nf <- matrix(rep(sf,each=nrow(object)),ncol=ncol(object))
    }
    # if the intercept is not finite, these rows
    # were all zero. here we put a small value instead
    intercept <- as.numeric(intercept)
    infiniteIntercept <- !is.finite(intercept)
    intercept[infiniteIntercept] <- -10
    normalizationFactors(object) <- nf * 2^intercept
    # we set the intercept, so replace the all zero
    # column with the rows which were all zero
    # in the previous dataset
    mcols(object)$allZero <- infiniteIntercept
  }

  # only continue on the rows with non-zero row sums
  objectNZ <- object[!mcols(object)$allZero,]
  stopifnot(all(!is.na(mcols(objectNZ)$dispFit)))
  
  # if a prior sigma squared not provided, estimate this
  # by the matching upper quantiles of the
  # log2 counts plus a pseudocount
  if (missing(betaPriorVar)) {
    logCounts <- log2(counts(objectNZ,normalized=TRUE) + 0.5)
    logFoldChangeMatrix <- logCounts - log2(mcols(objectNZ)$baseMean + 0.5)
    logFoldChangeVector <- as.numeric(logFoldChangeMatrix)
    varlogk <- 1/mcols(objectNZ)$baseMean + mcols(objectNZ)$dispFit
    weights <- 1/varlogk   
    betaPriorVar <- matchWeightedUpperQuantileForVariance(logFoldChangeVector, rep(weights,ncol(objectNZ)))
  }
  stopifnot(length(betaPriorVar)==1)
  
  lambda <- 1/rep(betaPriorVar,ncol(modelMatrix))
  # except for intercept which we set to wide prior
  if ("Intercept" %in% modelMatrixNames) {
    lambda[which(modelMatrixNames == "Intercept")] <- 1e-6
  }
  
  fit <- fitNbinomGLMs(object=objectNZ, modelMatrix=modelMatrix,
                       lambda=lambda, renameCols=FALSE,
                       alpha_hat=mcols(objectNZ)$dispFit,
                       betaTol=1e-4, useOptim=FALSE,
                       useQR=TRUE)
  normalizedDataNZ <- t(modelMatrix %*% t(fit$betaMatrix))

  normalizedData <- buildMatrixWithZeroRows(normalizedDataNZ, mcols(object)$allZero)

  # add back in the intercept, if finite
  if (!missing(intercept)) {
    normalizedData <- normalizedData + ifelse(infiniteIntercept, 0, intercept)
  }
  colnames(normalizedData) <- colnames(object)
  attr(normalizedData,"betaPriorVar") <- betaPriorVar
  if ("Intercept" %in% modelMatrixNames) {
    fittedInterceptNZ <- fit$betaMatrix[,which(modelMatrixNames == "Intercept"),drop=FALSE]
    fittedIntercept <- buildMatrixWithNARows(fittedInterceptNZ, mcols(object)$allZero)
    fittedIntercept[is.na(fittedIntercept)] <- -Inf
    attr(normalizedData,"intercept") <- fittedIntercept
  }
  normalizedData
}

rlogDataFast <- function(object, intercept, betaPriorVar, B) {
  if (is.null(mcols(object)$dispFit)) {
    stop("first estimate dispersion")
  }
  if (is.null(mcols(object)$allZero) | is.null(mcols(object)$baseMean)) {
    object <- getBaseMeansAndVariances(object)
  }
  if (!missing(B)) {
    if (length(B) != nrow(object)) {
      stop("B should be as long as the number of rows of object")
    }
    if (!all(B >= 0 & B <= 1)) {
      stop("B should be defined between 0 and 1")
    }
  }
  if (missing(intercept)) {
    # set the intercept to log2 ( geometric mean of normalized counts )
    intercept <- rowMeans(log2(counts(object,normalized=TRUE) + 0.5))
    interceptOut <- intercept
    interceptOut[mcols(object)$allZero] <- -Inf
  } else {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
    infiniteIntercept <- !is.finite(intercept)
    mcols(object)$allZero <- infiniteIntercept
    interceptOut <- intercept
  }
  # only continue on the rows with non-zero row sums
  objectNZ <- object[!mcols(object)$allZero,]
  stopifnot(all(!is.na(mcols(objectNZ)$dispFit)))
  interceptNZ <- intercept[!mcols(object)$allZero]

  logCounts <- log2(counts(objectNZ,normalized=TRUE) + 0.5)
  logFoldChangeMatrix <- logCounts - interceptNZ
  
  if (missing(B)) {
    # if a prior sigma squared notq provided, estimate this
    # by the matching upper quantiles of the
    # log2 counts plus a pseudocount
    
    if (missing(betaPriorVar)) {
      betaPriorVar <- if (nrow(objectNZ) > 1000) {
        # subsample 1000 genes along the range of base mean
        # to speed up the weighted quantile matching      
        idx <- order(mcols(objectNZ)$baseMean)[round(seq(from=1,to=nrow(objectNZ),length=1000))]
        logFoldChangeVector <- as.numeric(logFoldChangeMatrix[idx,,drop=FALSE])
        varlogk <- 1/mcols(objectNZ)$baseMean[idx] + mcols(objectNZ)$dispFit[idx]
        weights <- 1/varlogk
        matchWeightedUpperQuantileForVariance(logFoldChangeVector, rep(weights,ncol(objectNZ)))        
      } else {
        logFoldChangeVector <- as.numeric(logFoldChangeMatrix)
        varlogk <- 1/mcols(objectNZ)$baseMean + mcols(objectNZ)$dispFit
        weights <- 1/varlogk
        matchWeightedUpperQuantileForVariance(logFoldChangeVector, rep(weights,ncol(objectNZ)))
      }
    }
    stopifnot(length(betaPriorVar)==1)
    if (!is.null(normalizationFactors(object))) {
      nf <- normalizationFactors(objectNZ)
    } else {
      sf <- sizeFactors(object)
      nf <- matrix(rep(sf,each=nrow(objectNZ)),ncol=ncol(objectNZ))
    }
    dispersion <- mcols(objectNZ)$dispFit

    delta <- .05
    bgrid <- seq(from=0, to=1, by=delta)

    # find a subset of rows which covers the dynamic range
    lbm <- log(mcols(objectNZ)$baseMean)
    idx <- if (nrow(objectNZ) > 1000) {
      order(lbm)[round(seq(from=1, to=nrow(objectNZ), length=1000))]
    } else {
      seq_len(nrow(objectNZ))
    }

    # evaluate over a grid of B (shrinkage amount)
    optimalFromGrid <- rlogGridWrapper(counts(objectNZ)[idx,], nf[idx,], logFoldChangeMatrix[idx,],
                                       dispersion[idx], interceptNZ[idx], bgrid, betaPriorVar)
    fit <- loess(optimalFromGrid ~ lbm, data=data.frame(optimalFromGrid,lbm=lbm[idx]),
                 control=loess.control(trace.hat="approximate"), span=.1)
    pred <- predict(fit, newdata=data.frame(lbm))
    optimalB <- pmax(0, pmin(1, pred))
     
  } else {
    Bout <- B
    optimalB <- B[!mcols(object)$allZero]
  }

  # shrink the betas/log fold changes according to the interpolated B's
  lfcShrink <- (1 - optimalB) * logFoldChangeMatrix
  normalizedDataNZ <- interceptNZ + lfcShrink
  normalizedData <- buildMatrixWithZeroRows(normalizedDataNZ, mcols(object)$allZero)
  colnames(normalizedData) <- colnames(object)
  if (!missing(betaPriorVar)) {
    attr(normalizedData,"betaPriorVar") <- betaPriorVar
  }
  attr(normalizedData,"intercept") <- interceptOut
  Bout <- buildNumericWithZeroRows(optimalB, mcols(object)$allZero)
  attr(normalizedData,"B") <- Bout
  normalizedData
}

# Evaluate the likelihood over a grid of B's
rlogGridWrapper <- function (ySEXP, nfSEXP, betaSEXP, alphaSEXP, interceptSEXP, bgridSEXP, betapriorvarSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(rlogGridWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to rlogGrid, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  Bvec <- rlogGrid(ySEXP=ySEXP, nfSEXP=nfSEXP, betaSEXP=betaSEXP, alphaSEXP=alphaSEXP, interceptSEXP=interceptSEXP,
                   bgridSEXP=bgridSEXP, betapriorvarSEXP=betapriorvarSEXP)$Bvec
  as.numeric(Bvec)
}

sparseTest <- function(x, p, t1, t2) {
  rs <- rowSums(x)
  rmx <- apply(x, 1, max)
  if (all(rs <= t1)) return(invisible())
  prop <- (rmx/rs)[rs > t1]
  total <- mean(prop > p)
  if (total > t2) warning("the rlog assumes that data is close to a negative binomial distribution, an assumption
which is sometimes not compatible with datasets where many genes have many zero counts
despite a few very large counts.
In this data, for ",round(total,3)*100,"% of genes with a sum of normalized counts above ",t1,", it was the case 
that a single sample's normalized count made up more than ",p*100,"% of the sum over all samples.
the threshold for this warning is ",t2*100,"% of genes. See plotSparsity(dds) for a visualization of this.
We recommend instead using the varianceStabilizingTransformation or shifted log (see vignette).")
}
