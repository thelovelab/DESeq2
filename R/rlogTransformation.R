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
#' The prior width is calculated as follows: coefficients are fit for a model
#' with a term for each sample and for the intercept. This would typically
#' result in an unidentifiable solution, so a very wide prior is used.
#' Then the prior variance is estimated by taking the mean of the
#' row-wise variance of the sample coefficients. A second and final GLM fit
#' is performed using this prior.
#' It is also possible to supply the variance of the prior.
#' See the vignette for an example of the use and a comparison with
#' \code{varianceStabilizingTransformation}
#'
#' @aliases rlogTransformation rlogData
#'
#' @param object a DESeqDataSet
#' @param blind logical, whether to blind the transformation to the experimental
#' design. blind=TRUE should be used for comparing samples in an manner unbiased by
#' prior information on samples, for example to perform sample QA (quality assurance).
#' blind=FALSE should be used for transforming data for downstream analysis,
#' where the full use of the design information should be made.
#' @param samplesVector a character vector or factor of the sample identifiers
#' @param priorSigmasq a single value, the variance of the prior on the sample betas,
#' which if missing is estimated from the rows which do not have any
#' zeros
#' @param rowVarQuantile the quantile of the row variances of log fold changes
#' which will be used to set the width of the prior
#' 
#' @return for \code{rlogTransformation},
#' a SummarizedExperiment with assay data elements equal to
#' \eqn{\log_2(q_{ij}) = X_{j.} \beta_i}{log2(q_ij) = X_j. * beta_i},
#' see formula at \code{\link{DESeq}}.
#' for \code{rlogData}, a \code{matrix} of the same dimension as the
#' count data, containing the transformed values.  
#'
#' @seealso \code{\link{plotPCA}}, \code{\link{varianceStabilizingTransformation}}
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(betaSd=1)
#' rld <- rlogTransformation(dds, blind=TRUE)
#' dists <- dist(t(assay(rld)))
#' plot(hclust(dists))
#'
#' @export
rlogTransformation <- function(object, blind=TRUE, samplesVector, priorSigmasq, rowVarQuantile=.9) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~ 1
    object <- estimateDispersions(object)
  }
  if (is.null(dispersions(object))) {
    object <- estimateDispersions(object)
  }
  SummarizedExperiment(
    assays = rlogData(object, samplesVector, priorSigmasq, rowVarQuantile),
    colData = colData(object),
    rowData = rowData(object),
    exptData = exptData(object))
}

#' @rdname rlogTransformation
#' @export
rlogData <- function(object, samplesVector, priorSigmasq, rowVarQuantile=.9) {
  if (is.null(dispersions(object))) {
    stop("first estimate dispersion with a design of formula(~ 1)")
  }
  if (missing(samplesVector)) {
    samplesVector <- as.character(seq_len(ncol(object)))
  }
  
  # make a design matrix with a term for every sample
  # this would typically produce unidentifiable solution
  # for the GLM, but we add priors for all terms except
  # the intercept
  samplesVector <- factor(samplesVector,levels=unique(samplesVector))
  samples <- factor(c("null_level",as.character(samplesVector)),
                   levels=c("null_level",levels(samplesVector)))
  modelMatrix <- model.matrix(~samples)[-1,]
  modelMatrixNames <- colnames(modelMatrix)
  
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]

  # if a prior sigma squared not provided, calculate it
  # from betas calculated with a wide prior
  if (missing(priorSigmasq)) {
    lambda <- rep(1e-4, ncol(modelMatrix))
    if ("(Intercept)" %in% modelMatrixNames) {
      lambda[which(modelMatrixNames == "(Intercept)")] <- 1e-6
    }    
    fit <- fitNbinomGLMs(objectNZ,modelMatrix=modelMatrix,lambda=lambda,renameCols=FALSE)
    # use rows which have no zeros
    useNoZeros <- apply(counts(objectNZ),1,function(x) all(x > 0))
    if (sum(useNoZeros) == 0) {
      stop("no rows found without zeros")
    } 
    # calculate priors on sample betas
    # take row means of squares of sample betas
    betaRowMeanSquared <- rowMeans(fit$betaMatrix[,-which(fit$modelMatrixNames == "Intercept")]^2)
    priorSigmasq <- quantile(betaRowMeanSquared[useNoZeros], rowVarQuantile)
  }
  
  lambda <- 1/rep(priorSigmasq,ncol(modelMatrix))
  # except for intercept which we set to wide prior
  if ("Intercept" %in% fit$modelMatrixNames) {
    lambda[which(fit$modelMatrixNames == "Intercept")] <- 1e-6
  }
  fit <- fitNbinomGLMs(objectNZ,modelMatrix=modelMatrix,lambda=lambda,renameCols=FALSE)
  normalizedDataNZ <- t(modelMatrix %*% t(fit$betaMatrix))
  normalizedData <- buildMatrixWithNARows(normalizedDataNZ, mcols(object)$allZero)
  colnames(normalizedData) <- colnames(object)
  return(normalizedData)
}
