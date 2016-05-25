#' Apply a variance stabilizing transformation (VST) to the count data
#'
#' This function calculates a variance stabilizing transformation (VST) from the
#' fitted dispersion-mean relation(s) and then transforms the count data (normalized
#' by division by the size factors or normalization factors), yielding a matrix
#' of values which are now approximately homoskedastic (having constant variance along the range
#' of mean values). The transformation also normalizes with respect to library size.
#' The \code{\link{rlog}} is less sensitive
#' to size factors, which can be an issue when size factors vary widely.
#' These transformations are useful when checking for outliers or as input for
#' machine learning techniques such as clustering or linear discriminant analysis.
#' 
#' @aliases varianceStabilizingTransformation getVarianceStabilizedData
#' 
#' @param object a DESeqDataSet or matrix of counts
#' @param blind logical, whether to blind the transformation to the experimental
#' design. blind=TRUE should be used for comparing samples in an manner unbiased by
#' prior information on samples, for example to perform sample QA (quality assurance).
#' blind=FALSE should be used for transforming data for downstream analysis,
#' where the full use of the design information should be made.
#' blind=FALSE will skip re-estimation of the dispersion trend, if this has already been calculated.
#' If many of genes have large differences in counts due to
#' the experimental design, it is important to set blind=FALSE for downstream
#' analysis.
#' @param fitType in case dispersions have not yet been estimated for \code{object},
#' this parameter is passed on to \code{\link{estimateDispersions}} (options described there).
#'
#' @details For each sample (i.e., column of \code{counts(dds)}), the full variance function
#' is calculated from the raw variance (by scaling according to the size factor and adding 
#' the shot noise). We recommend a blind estimation of the variance function, i.e.,
#' one ignoring conditions. This is performed by default, and can be modified using the
#' 'blind' argument.
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
#' A typical workflow is shown in Section \emph{Variance stabilizing transformation}
#' in the package vignette.
#'
#' If \code{\link{estimateDispersions}} was called with:
#'
#' \code{fitType="parametric"},
#' a closed-form expression for the variance stabilizing
#' transformation is used on the normalized
#' count data. The expression can be found in the file \file{vst.pdf}
#' which is distributed with the vignette.
#'
#' \code{fitType="local"},
#' the reciprocal of the square root of the variance of the normalized counts, as derived
#' from the dispersion fit, is then numerically
#' integrated, and the integral (approximated by a spline function) is evaluated for each
#' count value in the column, yielding a transformed value. 
#'
#' \code{fitType="mean"}, a VST is applied for Negative Binomial distributed counts, 'k',
#' with a fixed dispersion, 'a': ( 2 asinh(sqrt(a k)) - log(a) - log(4) )/log(2).
#' 
#' In all cases, the transformation is scaled such that for large
#' counts, it becomes asymptotically (for large values) equal to the
#' logarithm to base 2 of normalized counts.
#'
#' The variance stabilizing transformation from a previous dataset
#' can be frozen and reapplied to new samples. See the 'Data quality assessment'
#' section of the vignette for strategies to see if new samples are
#' sufficiently similar to previous datasets. 
#' The frozen VST is accomplished by saving the dispersion function
#' accessible with \code{\link{dispersionFunction}}, assigning this
#' to the \code{DESeqDataSet} with the new samples, and running
#' varianceStabilizingTransformation with 'blind' set to FALSE
#' (see example below).
#' Then the dispersion function from the previous dataset will be used
#' to transform the new sample(s).
#'  
#' Limitations: In order to preserve normalization, the same
#' transformation has to be used for all samples. This results in the
#' variance stabilizition to be only approximate. The more the size
#' factors differ, the more residual dependence of the variance on the
#' mean will be found in the transformed data. \code{\link{rlog}} is a
#' transformation which can perform better in these cases.
#' As shown in the vignette, the function \code{meanSdPlot}
#' from the package \pkg{vsn} can be used to see whether this is a problem.
#'
#' @return \code{varianceStabilizingTransformation} returns a
#' \code{\link{DESeqTransform}} if a \code{DESeqDataSet} was provided,
#' or returns a a matrix if a count matrix was provided.
#' Note that for \code{\link{DESeqTransform}} output, the matrix of
#' transformed values is stored in \code{assay(vsd)}.
#' \code{getVarianceStabilizedData} also returns a matrix.
#'
#' @references
#'
#' Reference for the variance stabilizing transformation for counts with a dispersion trend:
#' 
#' Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 2010, 11:106. \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#' 
#' @author Simon Anders
#'
#' @seealso \code{\link{plotPCA}}, \code{\link{rlog}}, \code{\link{normTransform}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=6)
#' vsd <- varianceStabilizingTransformation(dds)
#' dists <- dist(t(assay(vsd)))
#' plot(hclust(dists))
#'
#' # learn the dispersion function of a dataset
#' design(dds) <- ~ 1
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#'
#' # use the previous dispersion function for a new sample
#' ddsNew <- makeExampleDESeqDataSet(m=1)
#' ddsNew <- estimateSizeFactors(ddsNew)
#' dispersionFunction(ddsNew) <- dispersionFunction(dds)
#' vsdNew <- varianceStabilizingTransformation(ddsNew, blind=FALSE)
#' 
#' @export
varianceStabilizingTransformation <- function (object, blind=TRUE, fitType="parametric") {
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names=colnames(object)), ~1)
  } else {
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (blind) {
    design(object) <- ~ 1
  }
  if (blind | is.null(attr(dispersionFunction(object),"fitType"))) {
    object <- estimateDispersionsGeneEst(object, quiet=TRUE)
    object <- estimateDispersionsFit(object, quiet=TRUE, fitType)
  }
  vsd <- getVarianceStabilizedData(object)
  if (matrixIn) {
    return(vsd)
  }
  se <- SummarizedExperiment(
    assays = vsd,
    colData = colData(object),
    rowRanges = rowRanges(object),
    metadata = metadata(object))
  DESeqTransform(se)
}

#' @rdname varianceStabilizingTransformation
#' @export
getVarianceStabilizedData <- function(object) {
  if (is.null(attr(dispersionFunction(object),"fitType"))) {
    stop("call estimateDispersions before calling getVarianceStabilizedData")
  }
  ncounts <- counts(object, normalized=TRUE)
  if( attr( dispersionFunction(object), "fitType" ) == "parametric" ) {
    coefs <- attr( dispersionFunction(object), "coefficients" )
    vst <- function( q ) {
      log( (1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * q + 2 * sqrt( coefs["asymptDisp"] * q * ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * q ) ) ) / ( 4 * coefs["asymptDisp"] ) ) / log(2)
    }
    return(vst(ncounts))
  } else if ( attr( dispersionFunction(object), "fitType" ) == "local" ) {
    # non-parametric fit -> numerical integration
    if (is.null(sizeFactors(object))) {
      stopifnot(!is.null(normalizationFactors(object)))
      # approximate size factors from columns of NF
      sf <- exp(colMeans(log(normalizationFactors(object))))
    } else {
      sf <- sizeFactors(object)
    }
    xg <- sinh( seq( asinh(0), asinh(max(ncounts)), length.out=1000 ) )[-1]
    xim <- mean( 1/sf )
    baseVarsAtGrid <- dispersionFunction(object)( xg ) * xg^2 + xim * xg
    integrand <- 1 / sqrt( baseVarsAtGrid )
    splf <- splinefun(
      asinh( ( xg[-1] + xg[-length(xg)] )/2 ),
      cumsum(
        ( xg[-1] - xg[-length(xg)] ) *
        ( integrand[-1] + integrand[-length(integrand)] )/2 ) )
    h1 <- quantile( rowMeans(ncounts), .95 )
    h2 <- quantile( rowMeans(ncounts), .999 )
    eta <- ( log2(h2) - log2(h1) ) / ( splf(asinh(h2)) - splf(asinh(h1)) )
    xi <- log2(h1) - eta * splf(asinh(h1))
    tc <- sapply( colnames(counts(object)), function(clm) {
      eta * splf( asinh( ncounts[,clm] ) ) + xi
    })
    rownames( tc ) <- rownames( counts(object) )
    return(tc)
  } else if ( attr( dispersionFunction(object), "fitType" ) == "mean" ) {
    alpha <- attr( dispersionFunction(object), "mean" )
    # the following stablizes NB counts with fixed dispersion alpha
    # and converges to log2(q) as q => infinity
    vst <- function(q) ( 2 * asinh(sqrt(alpha * q)) - log(alpha) - log(4) ) / log(2)
    return(vst(ncounts))
  } else {
    stop( "fitType is not parametric, local or mean" )
  }
}

#' Quickly estimate dispersion trend and apply a variance stabilizing transformation
#'
#' This is a wrapper for the \code{\link{varianceStabilizingTransformation}} (VST)
#' that provides much faster estimation of the dispersion trend used to determine
#' the formula for the VST. The speed-up is accomplished by
#' subsetting to a smaller number of genes in order to estimate this dispersion trend.
#' The subset of genes is chosen deterministically, to span the range
#' of genes' mean normalized count.
#' This wrapper for the VST is not blind to the experimental design:
#' the sample covariate information is used to estimate the global trend
#' of genes' dispersion values over the genes' mean normalized count.
#' It can be made strictly blind to experimental design by first
#' assigning a \code{\link{design}} of \code{~1} before running this function,
#' or by avoiding subsetting and using \code{\link{varianceStabilizingTransformation}}.
#' 
#' @param object a DESeqDataSet or a matrix of counts
#' @param blind logical, whether to blind the transformation to the experimental
#' design (see \code{\link{varianceStabilizingTransformation}})
#' @param nsub the number of genes to subset to (default 1000)
#' @param fitType for estimation of dispersions: this parameter
#' is passed on to \code{\link{estimateDispersions}} (options described there)
#'
#' @return a DESeqTranform object or a matrix of transformed, normalized counts
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(n=20000, m=20)
#' vsd <- vst(dds)
#'
#' @export
vst <- function(object, blind=TRUE, nsub=1000, fitType="parametric") {
  if (nrow(object) < nsub) {
    stop("less than 'nsub' rows,
  it is recommended to use varianceStabilizingTransformation directly")
  }
  if (is.null(colnames(object))) {
    colnames(object) <- seq_len(ncol(object))
  }
  if (is.matrix(object)) {
    matrixIn <- TRUE
    object <- DESeqDataSetFromMatrix(object, DataFrame(row.names=colnames(object)), ~ 1)
  } else {
    if (blind) {
      design(object) <- ~ 1
    }
    matrixIn <- FALSE
  }
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  baseMean <- rowMeans(counts(object, normalized=TRUE))
  if (sum(baseMean > 5) < nsub) {
    stop("less than 'nsub' rows with mean normalized count > 5, 
  it is recommended to use varianceStabilizingTransformation directly")
  }

  # subset to a specified number of genes with mean normalized count > 5
  object.sub <- object[baseMean > 5,]
  baseMean <- baseMean[baseMean > 5]
  o <- order(baseMean)
  idx <- o[round(seq(from=1, to=length(o), length=nsub))]
  object.sub <- object.sub[idx,]

  # estimate dispersion trend
  object.sub <- estimateDispersionsGeneEst(object.sub, quiet=TRUE)
  object.sub <- estimateDispersionsFit(object.sub, fitType=fitType, quiet=TRUE)

  # assign to the full object
  suppressMessages({dispersionFunction(object) <- dispersionFunction(object.sub)})

  # calculate and apply the VST
  vsd <- varianceStabilizingTransformation(object, blind=FALSE)
  if (matrixIn) {
    return(assay(vsd))
  } else {
    return(vsd)
  }
}
