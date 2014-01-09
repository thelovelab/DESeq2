#' Apply a 'regularized log' transformation
#'
#' This function uses Tikhonov/ridge regularization, as in \code{\link{nbinomWaldTest}},
#' to transform the data to the log2 scale in a way 
#' which minimizes differences between samples for rows with small counts.
#' The transformation produces a similar variance stabilizing effect as
#' \code{\link{varianceStabilizingTransformation}},
#' though \code{rlogTransformation} is more robust in the
#' case when the size factors vary widely.
#' The transformation is useful when checking for outliers
#' or as input for machine learning techniques
#' such as clustering or linear discriminant analysis.
#'
#' The 'regularization' referred to here corresponds to the maximum a posteriori
#' solution to the GLM with a prior on the coefficients for each sample.
#' The fitted dispersions are used rather than the MAP dispersions
#' (so similar to the \code{\link{varianceStabilizingTransformation}})
#' as the blind dispersion estimation would otherwise shrink
#' large, true log fold changes.
#' The prior variance is calculated as follows: 
#' a matrix is constructed of the logarithm of the counts plus a pseudocount of 0.5,
#' the log of the row means is then subtracted, leaving an estimate of
#' the log fold changes per sample over the fitted value using only an intercept.
#' The prior variance is then calculated as with
#' \code{\link{nbinomWaldTest}}, by matching the upper quantiles of the observed 
#' log fold change estimates with an upper quantile of the normal distribution.
#' A second and final GLM fit is then performed using this prior.
#' It is also possible to supply the variance of the prior.
#' See the vignette for an example of the use and a comparison with
#' \code{varianceStabilizingTransformation}.
#'
#' The transformed values, rlog(K), are equal to
#' \eqn{rlog(K_{ij}) = \log_2(q_{ij}) = x_{j.} \beta_i}{rlog(K_ij) = log2(q_ij) = x_j. * beta_i},
#' with formula terms defined in \code{\link{DESeq}}.
#'
#' The parameters of the rlog transformation from a previous dataset
#' can be "frozen" and reapplied to new samples. See the "Data quality assessment"
#' section of the vignette for strategies to see if new samples are
#' sufficiently similar to previous datasets. 
#' The "freezing" is accomplished by saving the dispersion function,
#' beta prior variance and the intercept from a previous dataset,
#' and running rlogTransformation with 'blind' set to FALSE
#' (see example below).
#' 
#' @aliases rlogTransformation rlogData rlogDataFast
#'
#' @param object a DESeqDataSet
#' @param blind logical, whether to blind the transformation to the experimental
#' design. blind=TRUE should be used for comparing samples in an manner unbiased by
#' prior information on samples, for example to perform sample QA (quality assurance).
#' blind=FALSE should be used for transforming data for downstream analysis,
#' where the full use of the design information should be made.
#' @param fast if TRUE, an approximation to the rlog transformation is made,
#' wherein an optimal amount of shrinkage of sample-wise log fold changes
#' is calculated for genes divided into 12 batched based on mean counts.
#' then shrinkage is performed by interpolating between these optimal values. 
#' optimal is defined by maximizing the same posterior as in the standard
#' rlog transformation.
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
#' 
#' @return for \code{rlogTransformation}, a \code{\link{SummarizedExperiment}}.
#' The matrix of transformed values are accessed by \code{assay(rld)}.
#' for \code{rlogData}, a \code{matrix} of the same dimension as the
#' count data, containing the transformed values. To avoid returning
#' matrices with NA values where there were zeros for all rows of
#' the unnormalized counts, rlogTransformation returns instead all
#' zeros (essentially adding a pseudocount of one, only to those rows 
#' in which all samples have zero).
#'
#' @seealso \code{\link{plotPCA}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=6,betaSD=1)
#' rld <- rlogTransformation(dds, blind=TRUE)
#' dists <- dist(t(assay(rld)))
#' plot(hclust(dists))
#'
#' # run the rlog transformation on one dataset
#' design(dds) <- ~ 1
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' rld <- rlogTransformation(dds, blind=FALSE)
#'
#' # apply the parameters to a new sample
#' 
#' ddsNew <- makeExampleDESeqDataSet(m=1)
#' mcols(ddsNew)$dispFit <- mcols(dds)$dispFit
#' betaPriorVar <- attr(rld,"betaPriorVar")
#' intercept <- mcols(rld)$rlogIntercept
#' rldNew <- rlogTransformation(ddsNew, blind=FALSE,
#'                              intercept=intercept,
#'                              betaPriorVar=betaPriorVar)
#'                            
#' 
#' @export
rlogTransformation <- function(object, blind=TRUE, fast=FALSE,
                               intercept, betaPriorVar, B) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~ 1
    object <- estimateDispersionsGeneEst(object)
    object <- estimateDispersionsFit(object)
  }
  if (is.null(mcols(object)$dispFit)) {
    object <- estimateDispersionsGeneEst(object)
    object <- estimateDispersionsFit(object)
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
    # if missing betaPriorVariance, just set a value to allow
    # the call below, it will be ignored within rlogDataFast()
    if (missing(betaPriorVar)) {
      betaPriorVar <- 1
    }
  }
  if (fast) {
    rld <- rlogDataFast(object, intercept, betaPriorVar, B)
  } else {
    rld <- rlogData(object, intercept, betaPriorVar)
  }
  se <- SummarizedExperiment(
           assays = rld,
           colData = colData(object),
           rowData = rowData(object),
           exptData = exptData(object))
  attr(se,"betaPriorVar") <- attr(rld, "betaPriorVar")
  if (!is.null(attr(rld,"intercept"))) {
    mcols(se)$rlogIntercept <- attr(rld,"intercept")
  }
  if (!is.null(attr(rld,"B"))) {
    attr(se,"B") <- attr(rld,"B")
  }
  se
}

#' @rdname rlogTransformation
#' @export
rlogData <- function(object, intercept, betaPriorVar) {
  if (is.null(mcols(object)$dispFit)) {
    stop("first estimate dispersion with a design of formula(~ 1)")
  }
  samplesVector <- as.character(seq_len(ncol(object)))
  if (!missing(intercept)) {
    if (length(intercept) != nrow(object)) {
      stop("intercept should be as long as the number of rows of object")
    }
  }
  if (!"allZero" %in% names(mcols(object))) {
    mcols(object)$allZero <- rowSums(counts(object)) == 0
  }
  if (!"baseMean" %in% names(mcols(object))) {
    mcols(object)$baseMean <- rowMeans(counts(object,normalized=TRUE))
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
    baseMean <- rowMeans(counts(objectNZ,normalized=TRUE))
    logFoldChangeMatrix <- logCounts - log2(baseMean + 0.5)
    logFoldChangeVector <- as.numeric(logFoldChangeMatrix)
    betaPriorVar <- matchUpperQuantileForVariance(logFoldChangeVector)
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



#' @rdname rlogTransformation
#' @export
rlogDataFast <- function(object, intercept, betaPriorVar, B) {
  if (is.null(mcols(object)$dispFit)) {
    stop("first estimate dispersion with a design of formula(~ 1)")
  }
  if (!"allZero" %in% names(mcols(object))) {
    mcols(object)$allZero <- rowSums(counts(object)) == 0
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
      logFoldChangeVector <- as.numeric(logFoldChangeMatrix)
      betaPriorVar <- matchUpperQuantileForVariance(logFoldChangeVector)
    }
    stopifnot(length(betaPriorVar)==1)
    if (!is.null(normalizationFactors(object))) { 
      nf <- normalizationFactors(objectNZ)
    } else {
      sf <- sizeFactors(object)
      nf <- matrix(rep(sf,each=nrow(objectNZ)),ncol=ncol(objectNZ))
    }
    dispersion <- mcols(objectNZ)$dispFit
    
    # the posterior / penalized log likelihood times -1
    objectiveFn <- function(B, idx) {
      lfcShrink <- (1 - B) * logFoldChangeMatrix
      logLike <- sum(dnbinom(counts(objectNZ)[idx,], mu = (nf * 2^(interceptNZ + lfcShrink))[idx,],
                             size = 1/dispersion[idx], log=TRUE))
      logPrior <- sum(dnorm(lfcShrink[idx,], 0, sqrt(betaPriorVar), log=TRUE))
      -1 * (logLike + logPrior)
    }

    nq <- 12
    quantiles <- quantile(interceptNZ, 0:nq/nq)
    quantiles[1] <- quantiles[1] - 1
    quantiles[11] <- quantiles[11] + 1
    cutByQuantiles <- as.integer(cut(interceptNZ, quantiles))
    
    Bs <- sapply(1:nq, function(i) {
      idx <- cutByQuantiles == i
      optimize(objectiveFn, interval=c(0,1), idx=idx)$minimum
    })

    # interpolate between the optimal B's based on the log2(baseMean)
    xs <- (quantiles[-1] + quantiles[-(nq + 1)])/2
    logit <- function(p) log(p / (1 - p))
    sigmoid <- function(x) 1 / (1 + exp(-x))
    logitBs <- logit(Bs)
    lofit <- loess(logitBs ~ xs)
    pred <- predict(lofit, newdata=data.frame(xs=interceptNZ))
    pred[interceptNZ <= xs[1]] <- logitBs[1]
    pred[interceptNZ >= xs[nq]] <- logitBs[nq]
    optimalB <- sigmoid(pred)
  } else {
    Bout <- B
    optimalB <- B[!mcols(object)$allZero]
  }
  
  # shrink the betas/log fold changes according to the interpolated B's
  lfcShrink <- (1 - optimalB) * logFoldChangeMatrix
  normalizedDataNZ <- interceptNZ + lfcShrink
  normalizedData <- buildMatrixWithZeroRows(normalizedDataNZ, mcols(object)$allZero)
  colnames(normalizedData) <- colnames(object)
  attr(normalizedData,"betaPriorVar") <- betaPriorVar
  attr(normalizedData,"intercept") <- interceptOut
  Bout <- buildNumericWithZeroRows(optimalB, mcols(object)$allZero)
  attr(normalizedData,"B") <- Bout
  normalizedData
}
