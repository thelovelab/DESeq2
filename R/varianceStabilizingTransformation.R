#' Apply a variance stabilizing transformation (VST) to the count data
#'
#' This function calculates a variance stabilizing transformation (VST) from the
#' fitted dispersion-mean relation(s) and then transforms the count data (normalized
#' by division by the size factors or normalization factors), yielding a matrix
#' of values which are now approximately homoskedastic (having constant variance along the range
#' of mean values). The \code{\link{rlogTransformation}} is less sensitive
#' to size factors, which can be an issue when size factors vary widely.
#' This transformation is useful when checking for outliers or as input for
#' machine learning techniques such as clustering or linear discriminant analysis.
#'
#' @aliases varianceStabilizingTransformation getVarianceStabilizedData
#' 
#' @param object a DESeqDataSet, with \code{design(object) <- formula(~ 1)}
#' and size factors (or normalization factors) and dispersions estimated
#' using local or parametric \code{fitType}.
#'
#' @details For each sample (i.e., column of \code{counts(dds)}), the full variance function
#'   is calculated from the raw variance (by scaling according to the size factor and adding 
#'   the shot noise). The function requires a blind estimate of the variance function, i.e.,
#'   one ignoring conditions. This is achieved by setting the design formula to \code{~ 1} and
#'   then calling \code{\link{estimateDispersions}}.
#'
#'   A typical workflow is shown in Section \emph{Variance stabilizing transformation} in the package vignette.
#'
#'   If \code{\link{estimateDispersions}} was called with \code{fitType="parametric"},
#'   a closed-form expression for the variance stabilizing transformation is used on the normalized
#'   count data. The expression can be found in the file \file{vst.pdf} which is distributed with the vignette.
#'
#'   If \code{\link{estimateDispersions}} was called with \code{fitType="local"},
#'   the reciprocal of the square root of the variance of the normalized counts, as derived
#'   from the dispersion fit, is then numerically
#'   integrated, and the integral (approximated by a spline function) is evaluated for each
#'   count value in the column, yielding a transformed value. 
#'   
#'   In both cases, the transformation is scaled such that for large
#'   counts, it becomes asymptotically (for large values) equal to the
#'   logarithm to base 2.
#'
#'   Limitations: In order to preserve normalization, the same
#'   transformation has to be used for all samples. This results in the
#'   variance stabilizition to be only approximate. The more the size
#'   factors differ, the more residual dependence of the variance on the
#'   mean you will find in the transformed data. As shown in the vignette, you can use the function
#'   \code{meanSdPlot} from the package \pkg{vsn} to see whether this
#'   is a problem for your data.
#'
#' @return for \code{varianceStabilizingTransformation}, a \code{SummarizedExperiment}.
#' for \code{getVarianceStabilizedData}, a \code{matrix} of the same dimension as the
#' count data, containing the transformed values.  
#' 
#' @author Simon Anders
#'
#' @seealso \code{\link{plotPCA}}, \code{\link{rlogTransformation}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' design(dds) <- formula(~ 1)
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' vsd <- varianceStabilizingTransformation(dds)
#' par(mfrow=c(1,2))
#' plot(rank(rowMeans(counts(dds))), genefilter::rowVars(log2(counts(dds)+1)), main="log2(x+1) transform")
#' plot(rank(rowMeans(assay(vsd))), genefilter::rowVars(assay(vsd)), main="VST")
#' 
#' @export
varianceStabilizingTransformation <- function (object) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    object <- estimateSizeFactors(object)
  }
  if (is.null(dispersions(object))) {
    object <- estimateDispersions(object)
  }
  SummarizedExperiment(
    assays = getVarianceStabilizedData(object),
    colData = colData(object),
    rowData = rowData(object),
    exptData = exptData(object))
}

#' @rdname varianceStabilizingTransformation
#' @export
getVarianceStabilizedData <- function(object) {
  if (is.null(attr(dispersionFunction(object),"fitType"))) {
    stop("call estimateDispersions before calling getVarianceStabilizedData")
  }
  ncounts <- counts(object,normalized=TRUE)
  if( attr( dispersionFunction(object), "fitType" ) == "parametric" ) {
    coefs <- attr( dispersionFunction(object), "coefficients" )
    vst <- function( q ) {
      log( (1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * q + 2 * sqrt( coefs["asymptDisp"] * q * ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * q ) ) ) / ( 4 * coefs["asymptDisp"] ) ) / log(2)
    }
    return(vst(ncounts))
  } else {
    # non-parametric fit -> numerical integration
    if (is.null(sizeFactors(object))) {
      stop("call estimateSizeFactors before calling getVarianceStabilizedData if using local dispersion fit")
    }
    xg <- sinh( seq( asinh(0), asinh(max(ncounts)), length.out=1000 ) )[-1]
    xim <- mean( 1/sizeFactors(object) )
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
  }
}
