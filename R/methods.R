counts.DESeqDataSet <- function(object, normalized=FALSE, replaced=FALSE) {
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  if (replaced) {
    if ("replaceCounts" %in% assayNames(object)) {
      cnts <- assays(object)[["replaceCounts"]]
    } else {
      warning("there are no assays named 'replaceCounts', using original.
calling DESeq() will replace outliers if they are detected and store this assay.")
      cnts <- assays(object)[["counts"]]
    }
  } else {
    cnts <- assays(object)[["counts"]]
  }
  if (!normalized) {
    return(cnts)
  } else {
    if (!is.null(normalizationFactors(object))) {
      return( cnts / normalizationFactors(object) )
    } else if (is.null(sizeFactors(object)) || any(is.na(sizeFactors(object)))) {
      stop("first calculate size factors, add normalizationFactors, or set normalized=FALSE")
    } else {
      return( t( t( cnts ) / sizeFactors(object) ) )
    }
  }
}

#' Accessors for the 'counts' slot of a DESeqDataSet object.
#' 
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample. 
#'
#' @rdname counts
#'
#' @param object a \code{DESeqDataSet} object.
#' @param normalized logical indicating whether or not to divide the counts by
#' the size factors or normalization factors before returning
#' (normalization factors always preempt size factors)
#' @param replaced after a \code{DESeq} call, this argument will return the counts
#' with outliers replaced instead of the original counts, and optionally \code{normalized}.
#' The replaced counts are stored by \code{DESeq} in \code{assays(object)[['replaceCounts']]}.
#' @param value an integer matrix
#' @author Simon Anders
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalizationFactors}}
#'
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet(m=4)
#' head(counts(dds))
#'
#' dds <- estimateSizeFactors(dds) # run this or DESeq() first
#' head(counts(dds, normalized=TRUE))
#'
#' @export
setMethod("counts", signature(object="DESeqDataSet"), counts.DESeqDataSet)

#' @rdname counts
#' @export
setReplaceMethod("counts", signature(object="DESeqDataSet", value="matrix"),
                 function( object, value ) {
                   assays(object)[["counts"]] <- value
                   validObject(object)
                   object
                 })   


design.DESeqDataSet <- function(object) object@design

#' Accessors for the 'design' slot of a DESeqDataSet object.
#' 
#' The design holds the R \code{formula} which expresses how the
#' counts depend on the variables in \code{colData}.
#' See \code{\link{DESeqDataSet}} for details.
#' 
#' @rdname design
#'
#' @param object a \code{DESeqDataSet} object
#' @param value a \code{formula} used for estimating dispersion
#' and fitting Negative Binomial GLMs
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' design(dds) <- formula(~ 1)
#'
#' @export
setMethod("design", signature(object="DESeqDataSet"), design.DESeqDataSet)

design.replace <- function( object, value ) {
  # Temporary hack for backward compatibility with "old"
  # DESeqDataSet objects. Remove once all serialized
  # DESeqDataSet objects around have been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  object@design <- value
  validObject(object)
  object
}

#' @rdname design
#' @export
setReplaceMethod("design", signature(object="DESeqDataSet", value="formula"), design.replace)

#' @rdname design
#' @export
setReplaceMethod("design", signature(object="DESeqDataSet", value="matrix"), design.replace)

dispersionFunction.DESeqDataSet <- function(object) object@dispersionFunction

#' Accessors for the 'dispersionFunction' slot of a DESeqDataSet object.
#'
#' The dispersion function is calculated by \code{\link{estimateDispersions}} and
#' used by \code{\link{varianceStabilizingTransformation}}.  Parametric dispersion
#' fits store the coefficients of the fit as attributes in this slot.
#'
#' Setting this will also overwrite \code{mcols(object)$dispFit} and the estimate
#' the variance of dispersion residuals, see \code{estimateVar} below.
#'
#' @rdname dispersionFunction
#'
#' @param object a \code{DESeqDataSet} object.
#' @param value a \code{function}
#' @param ... additional arguments
#' 
#' @seealso \code{\link{estimateDispersions}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' dispersionFunction(dds)
#'
#' @export
setMethod("dispersionFunction", signature(object="DESeqDataSet"),
          dispersionFunction.DESeqDataSet)

dispFun.replace <- function(object, value) {
  # Temporary hack for backward compatibility with "old"
  # DESeqDataSet objects. Remove once all serialized
  # DESeqDataSet objects around have been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)

  # the following will add 'dispFit' to mcols(object)

  # first, check to see that we have 'baseMean' and 'allZero' columns
  if (is.null(mcols(object)$baseMean) | is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }
  # warning about existing 'dispFit' data will be removed
  if (!is.null(mcols(object)$dispFit)) {
    mcols(object) <- mcols(object)[,!names(mcols(object)) == "dispFit",drop=FALSE]
  }
  # now call the dispersionFunction on 'baseMean' to make 'dispFit'
  nonzeroIdx <- !mcols(object)$allZero
  dispFit <- value(mcols(object)$baseMean[nonzeroIdx])
  # if the function returns a single value, build the full vector
  if (length(dispFit) == 1) {
    dispFit <- rep(dispFit, sum(nonzeroIdx))
  }
  dispDataFrame <- buildDataFrameWithNARows(list(dispFit=dispFit),
                                            mcols(object)$allZero)
  mcols(dispDataFrame) <- DataFrame(type="intermediate",
                                    description="fitted values of dispersion")
  mcols(object) <- cbind(mcols(object), dispDataFrame)

  # estimate variance of log dispersion around the fit

  # need to estimate variance of log dispersion residuals
  minDisp <- 1e-8
  dispGeneEst <- mcols(object)$dispGeneEst[nonzeroIdx]
  aboveMinDisp <- dispGeneEst >= minDisp*100
  if (sum(aboveMinDisp,na.rm=TRUE) > 0) {
    dispResiduals <- log(dispGeneEst) - log(dispFit)
    varLogDispEsts <- mad(dispResiduals[aboveMinDisp],na.rm=TRUE)^2
    attr( value, "varLogDispEsts" ) <- varLogDispEsts
  } else {
    message("variance of dispersion residuals not estimated (necessary only for differential expression calling)")
  }

  # store the dispersion function
  object@dispersionFunction <- value
  validObject(object)
  object
}

#' @rdname dispersionFunction
#' @export
setReplaceMethod("dispersionFunction", signature(object="DESeqDataSet", value="function"), dispFun.replace)

dispersions.DESeqDataSet <- function(object) mcols(object)$dispersion

#' Accessor functions for the dispersion estimates in a DESeqDataSet
#' object.
#' 
#' The dispersions for each row of the DESeqDataSet.  Generally,
#' these are set by \code{\link{estimateDispersions}}.
#' 
#' @rdname dispersions
#'
#' @param object a \code{DESeqDataSet} object.
#' @param value the dispersions to use for the Negative Binomial modeling
#' @param ... additional arguments
#'
#' @author Simon Anders
#' @seealso \code{\link{estimateDispersions}}
#' 
#' @export
setMethod("dispersions", signature(object="DESeqDataSet"),
          dispersions.DESeqDataSet)

#' @rdname dispersions
#' @export
setReplaceMethod("dispersions", signature(object="DESeqDataSet", value="numeric"),
                 function(object, value) {
                   firstRowDataColumn <- ncol(mcols(object)) == 0
                   mcols(object)$dispersion <- value
                   if (firstRowDataColumn) {
                     mcols(mcols(object)) <- DataFrame(type="input",
                                                       description="final estimate of dispersion")
                   }
                   validObject( object )
                   object
                 })


sizeFactors.DESeqDataSet <- function(object) {
  if (!"sizeFactor" %in% names(colData(object))) return(NULL)
  sf <- object$sizeFactor
  names( sf ) <- colnames( object )
  sf
}

#' Accessor functions for the 'sizeFactors' information in a DESeqDataSet
#' object.
#' 
#' The sizeFactors vector assigns to each column of the count matrix a value, the
#' size factor, such that count values in the columns can be brought to a common
#' scale by dividing by the corresponding size factor (as performed by
#' \code{counts(dds, normalized=TRUE)}).
#' See \code{\link{DESeq}} for a description of the use of size factors.
#' If gene-specific normalization
#' is desired for each sample, use \code{\link{normalizationFactors}}.
#' 
#' @rdname sizeFactors
#'
#' @param object a \code{DESeqDataSet} object.
#' @param value a numeric vector, one size factor for each column in the count
#' data.
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactors}}
#'
#' @export
setMethod("sizeFactors", signature(object="DESeqDataSet"),
          sizeFactors.DESeqDataSet)

#' @rdname sizeFactors
#' @export
setReplaceMethod("sizeFactors", signature(object="DESeqDataSet", value="numeric"),
                 function( object, value ) {
                   stopifnot(all(!is.na(value)))
                   stopifnot(all(is.finite(value)))
                   stopifnot(all(value > 0))

                   # Temporary hack for backward compatibility with "old"
                   # DESeqDataSet objects. Remove once all serialized
                   # DESeqDataSet objects around have been updated.
                   if (!.hasSlot(object, "rowRanges"))
                     object <- updateObject(object)
                   # have to make sure to remove sizeFactor which might be
                   # coming from a previous CountDataSet
                   object$sizeFactor <- value
                   idx <- which(colnames(colData(object)) == "sizeFactor")
                   metaDataFrame <- DataFrame(type="intermediate",
                                              description="a scaling factor for columns")
                   mcols(colData(object))[idx,] <- metaDataFrame
                   validObject( object )
                   object
                 }) 

normalizationFactors.DESeqDataSet <- function(object) {
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  if (!"normalizationFactors" %in% assayNames(object)) return(NULL)
  assays(object)[["normalizationFactors"]]
}

#' Accessor functions for the normalization factors in a DESeqDataSet
#' object.
#'
#' Gene-specific normalization factors for each sample can be provided as a matrix,
#' which will preempt \code{\link{sizeFactors}}. In some experiments, counts for each
#' sample have varying dependence on covariates, e.g. on GC-content for sequencing
#' data run on different days, and in this case it makes sense to provide
#' gene-specific factors for each sample rather than a single size factor.
#'
#' Normalization factors alter the model of \code{\link{DESeq}} in the following way, for
#' counts \eqn{K_{ij}}{K_ij} and normalization factors \eqn{NF_{ij}}{NF_ij} for gene i and sample j:
#'
#' \deqn{ K_{ij} \sim \textrm{NB}( \mu_{ij}, \alpha_i) }{ K_ij ~ NB(mu_ij, alpha_i) }
#' \deqn{ \mu_{ij} = NF_{ij} q_{ij} }{ mu_ij = NF_ij q_ij }
#'
#' @note Normalization factors are on the scale of the counts (similar to \code{\link{sizeFactors}})
#' and unlike offsets, which are typically on the scale of the predictors (in this case, log counts).
#' Normalization factors should include library size normalization. They should have
#' row-wise geometric mean near 1, as is the case with size factors, such that the mean of normalized
#' counts is close to the mean of unnormalized counts. See example code below.
#'
#' @docType methods
#' @name normalizationFactors
#' @rdname normalizationFactors
#' @aliases normalizationFactors normalizationFactors,DESeqDataSet-method normalizationFactors<-,DESeqDataSet,matrix-method
#' @param object a \code{DESeqDataSet} object.
#' @param value the matrix of normalization factors
#' @param ... additional arguments
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(n=100, m=4)
#'
#' normFactors <- matrix(runif(nrow(dds)*ncol(dds),0.5,1.5),
#'                       ncol=ncol(dds),nrow=nrow(dds),
#'                       dimnames=list(1:nrow(dds),1:ncol(dds)))
#'
#' # the normalization factors matrix should not have 0's in it
#' # it should have geometric mean near 1 for each row
#' normFactors <- normFactors / exp(rowMeans(log(normFactors)))
#' normalizationFactors(dds) <- normFactors
#'
#' dds <- DESeq(dds)
#'
#' @export
setMethod("normalizationFactors", signature(object="DESeqDataSet"),
          normalizationFactors.DESeqDataSet)

#' @rdname normalizationFactors
#' @export
setReplaceMethod("normalizationFactors", signature(object="DESeqDataSet", value="matrix"),
                 function(object, value) {
                   stopifnot(all(!is.na(value)))
                   stopifnot(all(is.finite(value)))
                   stopifnot(all(value > 0))

                   # Temporary hack for backward compatibility with "old"
                   # DESeqDataSet objects. Remove once all serialized
                   # DESeqDataSet objects around have been updated.
                   if (!.hasSlot(object, "rowRanges"))
                     object <- updateObject(object)
                   # enforce same dimnames
                   dimnames(value) <- dimnames(object)
                   assays(object)[["normalizationFactors"]] <- value
                   validObject( object )
                   object
                 })

estimateSizeFactors.DESeqDataSet <- function(object, type=c("ratio","poscounts","iterate"),
                                             locfunc=stats::median,
                                             geoMeans, controlGenes, normMatrix, quiet=FALSE) {
  type <- match.arg(type, c("ratio","poscounts","iterate"))
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges")) {
    object <- updateObject(object)
  }
  object <- sanitizeColData(object)
  if (type == "iterate") {
    sizeFactors(object) <- estimateSizeFactorsIterate(object)
  } else {
    if (type == "poscounts") {
      geoMeanNZ <- function(x) {
        if (all(x == 0)) { 0 } else { exp( sum(log(x[x > 0])) / length(x) ) }
      }
      geoMeans <- apply(counts(object), 1, geoMeanNZ)
    }
    if ("avgTxLength" %in% assayNames(object)) {
      nm <- assays(object)[["avgTxLength"]]
      nm <- nm / exp(rowMeans(log(nm))) # divide out the geometric mean
      normalizationFactors(object) <- estimateNormFactors(counts(object),
                                                          normMatrix=nm,
                                                          locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
      if (!quiet) message("using 'avgTxLength' from assays(dds), correcting for library size")
    } else if (missing(normMatrix)) {
      sizeFactors(object) <- estimateSizeFactorsForMatrix(counts(object), locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
    } else {
      normalizationFactors(object) <- estimateNormFactors(counts(object),
                                                          normMatrix=normMatrix,
                                                          locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
      if (!quiet) message("using 'normMatrix', adding normalization factors which correct for library size")
    }
  }
  object
}

#' Estimate the size factors for a \code{\link{DESeqDataSet}}
#' 
#' This function estimates the size factors using the
#' "median ratio method" described by Equation 5 in Anders and Huber (2010).
#' The estimated size factors can be accessed using the accessor function \code{\link{sizeFactors}}.
#' Alternative library size estimators can also be supplied
#' using the assignment function \code{\link{sizeFactors<-}}.
#' 
#' Typically, the function is called with the idiom:
#'
#' \code{dds <- estimateSizeFactors(dds)}
#'
#' See \code{\link{DESeq}} for a description of the use of size factors in the GLM.
#' One should call this function after \code{\link{DESeqDataSet}}
#' unless size factors are manually specified with \code{\link{sizeFactors}}.
#' Alternatively, gene-specific normalization factors for each sample can be provided using
#' \code{\link{normalizationFactors}} which will always preempt \code{\link{sizeFactors}}
#' in calculations.
#'
#' Internally, the function calls \code{\link{estimateSizeFactorsForMatrix}}, 
#' which provides more details on the calculation.
#'
#' @rdname estimateSizeFactors
#' 
#' @param object a DESeqDataSet
#' @param type Method for estimation: either "ratio", "poscounts", or "iterate".
#' "ratio" uses the standard median ratio method introduced in DESeq. The size factor is the
#' median ratio of the sample over a "pseudosample": for each gene, the geometric mean
#' of all samples.
#' "poscounts" and "iterate" offer alternative estimators, which can be
#' used even when all genes contain a sample with a zero (a problem for the
#' default method, as the geometric mean becomes zero, and the ratio undefined).
#' The "poscounts" estimator deals with a gene with some zeros, by calculating a
#' modified geometric mean by taking the n-th root of the product of the non-zero counts.
#' This evolved out of use cases with Paul McMurdie's phyloseq package for metagenomic samples.
#' The "iterate" estimator iterates between estimating the dispersion with a design of ~1, and
#' finding a size factor vector by numerically optimizing the likelihood
#' of the ~1 model.
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used. However, especially for low counts, the
#' \code{\link[genefilter]{shorth}} function from the genefilter package may give better results.
#' @param geoMeans by default this is not provided and the
#' geometric means of the counts are calculated within the function.
#' A vector of geometric means from another count matrix can be provided
#' for a "frozen" size factor calculation. The size factors will be 
#' scaled to have a geometric mean of 1 when supplying \code{geoMeans}.
#' @param controlGenes optional, numeric or logical index vector specifying those genes to
#' use for size factor estimation (e.g. housekeeping or spike-in genes)
#' @param normMatrix optional, a matrix of normalization factors which do not yet
#' control for library size. Note that this argument should not be used (and
#' will be ignored) if the \code{dds} object was created using \code{tximport}.
#' In this case, the information in \code{assays(dds)[["avgTxLength"]]}
#' is automatically used to create appropriate normalization factors.
#' Providing \code{normMatrix} will estimate size factors on the
#' count matrix divided by \code{normMatrix} and store the product of the
#' size factors and \code{normMatrix} as \code{\link{normalizationFactors}}.
#' It is recommended to divide out the row-wise geometric mean of
#' \code{normMatrix} so the rows roughly are centered on 1.
#' @param quiet whether to print messages
#' 
#' @return The DESeqDataSet passed as parameters, with the size factors filled
#' in.
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactorsForMatrix}}
#'
#' @references
#'
#' Reference for the median ratio method:
#' 
#' Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data.
#' Genome Biology 2010, 11:106. \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#' 
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet(n=1000, m=4)
#' dds <- estimateSizeFactors(dds)
#' sizeFactors(dds)
#'
#' dds <- estimateSizeFactors(dds, controlGenes=1:200)
#'
#' m <- matrix(runif(1000 * 4, .5, 1.5), ncol=4)
#' dds <- estimateSizeFactors(dds, normMatrix=m)
#' normalizationFactors(dds)[1:3,]
#' 
#' geoMeans <- exp(rowMeans(log(counts(dds))))
#' dds <- estimateSizeFactors(dds,geoMeans=geoMeans)
#' sizeFactors(dds)
#'
#' @export
setMethod("estimateSizeFactors", signature(object="DESeqDataSet"),
          estimateSizeFactors.DESeqDataSet)

estimateDispersions.DESeqDataSet <- function(object, fitType=c("parametric","local","mean", "glmGamPoi"),
                                             maxit=100, useCR=TRUE,
                                             weightThreshold=1e-2,
                                             quiet=FALSE, modelMatrix=NULL,
                                             minmu=if (fitType=="glmGamPoi") 1e-6 else 0.5) {
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
  }
  # size factors could have slipped in to colData from a previous run
  if (!is.null(sizeFactors(object))) {
    if (!is.numeric(sizeFactors(object))) {
      stop("the sizeFactor column in colData is not numeric.
this column could have come in during colData import and should be removed.")
    }
    if (any(is.na(sizeFactors(object)))) {
      stop("the sizeFactor column in colData contains NA.
this column could have come in during colData import and should be removed.")
    }
  }
  if (all(rowSums(counts(object) == counts(object)[,1]) == ncol(object))) {
    stop("all genes have equal values for all samples. will not be able to perform differential analysis")
  }
  if (!is.null(dispersions(object))) {
    if (!quiet) message("found already estimated dispersions, replacing these")
    mcols(object) <- mcols(object)[,!(mcols(mcols(object))$type %in% c("intermediate","results")),drop=FALSE]
  }
  stopifnot(length(maxit)==1)
  fitType <- match.arg(fitType, choices=c("parametric","local","mean", "glmGamPoi"))
  dispersionEstimator <- if (fitType == "glmGamPoi") {
    if (!requireNamespace("glmGamPoi", quietly=TRUE)) {
      stop("type='glmGamPoi' requires installing the Bioconductor package 'glmGamPoi'")
    }
    "glmGamPoi"
  } else {
    "DESeq2"
  }
  
  
  checkForExperimentalReplicates(object, modelMatrix)
  
  if (!quiet) message("gene-wise dispersion estimates")
  object <- estimateDispersionsGeneEst(object,
                                       maxit=maxit,
                                       useCR=useCR,
                                       weightThreshold=weightThreshold,
                                       quiet=quiet,
                                       modelMatrix=modelMatrix,
                                       minmu=minmu,
                                       type = dispersionEstimator)
  if (!quiet) message("mean-dispersion relationship")
  object <- estimateDispersionsFit(object,
                                   fitType= fitType,
                                   quiet=quiet)
  if (!quiet) message("final dispersion estimates")
  object <- estimateDispersionsMAP(object,
                                   maxit=maxit,
                                   useCR=useCR,
                                   weightThreshold=weightThreshold,
                                   quiet=quiet,
                                   modelMatrix=modelMatrix,
                                   type = dispersionEstimator)
  
  return(object)
}

checkForExperimentalReplicates <- function(object, modelMatrix) {
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  
  noReps <- if (is.null(modelMatrix)) {
    mmtest <- getModelMatrix(object)
    nrow(mmtest) == ncol(mmtest)
  } else {
    nrow(modelMatrix) == ncol(modelMatrix)
  }
  if (noReps) {
    stop("

  The design matrix has the same number of samples and coefficients to fit,
  so estimation of dispersion is not possible. Treating samples
  as replicates was deprecated in v1.20 and no longer supported since v1.22.

")
  }
  TRUE
}

#' Estimate the dispersions for a DESeqDataSet
#' 
#' This function obtains dispersion estimates for Negative Binomial distributed data.
#'
#' Typically the function is called with the idiom:
#'
#' \code{dds <- estimateDispersions(dds)}
#'
#' The fitting proceeds as follows: for each gene, an estimate of the dispersion
#' is found which maximizes the Cox Reid-adjusted profile likelihood
#' (the methods of Cox Reid-adjusted profile likelihood maximization for
#' estimation of dispersion in RNA-Seq data were developed by McCarthy,
#' et al. (2012), first implemented in the edgeR package in 2010);
#' a trend line capturing the dispersion-mean relationship is fit to the maximum likelihood estimates;
#' a normal prior is determined for the log dispersion estimates centered
#' on the predicted value from the trended fit
#' with variance equal to the difference between the observed variance of the
#' log dispersion estimates and the expected sampling variance;
#' finally maximum a posteriori dispersion estimates are returned.
#' This final dispersion parameter is used in subsequent tests.
#' The final dispersion estimates can be accessed from an object using \code{\link{dispersions}}.
#' The fitted dispersion-mean relationship is also used in
#' \code{\link{varianceStabilizingTransformation}}.
#' All of the intermediate values (gene-wise dispersion estimates, fitted dispersion
#' estimates from the trended fit, etc.) are stored in \code{mcols(dds)}, with
#' information about these columns in \code{mcols(mcols(dds))}.
#'
#' The log normal prior on the dispersion parameter has been proposed
#' by Wu, et al. (2012) and is also implemented in the DSS package.
#'
#' In DESeq2, the dispersion estimation procedure described above replaces the
#' different methods of dispersion from the previous version of the DESeq package.
#' 
#' Since version 1.29, DESeq2 can call the glmGamPoi package, which can speed up the inference
#' and is optimized for fitting many samles with very small counts (for example single cell 
#' RNA-seq data). To call functions from the glmGamPoi package, make sure that it is installed
#' and set \code{fitType = "glmGamPoi"}. In addition, to the gene estimates, the trend and the MAP,
#' the glmGamPoi package calculates the corresponding quasi-likelihood estimates. Those can be
#' used with the \code{nbinomLRT()} test to get more precise p-value estimates.
#'
#' The lower-level functions called by \code{estimateDispersions} are:
#' \code{\link{estimateDispersionsGeneEst}},
#' \code{\link{estimateDispersionsFit}}, and
#' \code{\link{estimateDispersionsMAP}}.
#' 
#' @rdname estimateDispersions
#'
#' @param object a DESeqDataSet
#' @param fitType either "parametric", "local", "mean", or "glmGamPoi"
#' for the type of fitting of dispersions to the mean intensity.
#' \itemize{
#'   \item parametric - fit a dispersion-mean relation of the form:
#'     \deqn{dispersion = asymptDisp + extraPois / mean}
#'     via a robust gamma-family GLM. The coefficients \code{asymptDisp} and \code{extraPois}
#'     are given in the attribute \code{coefficients} of the \code{\link{dispersionFunction}}
#'     of the object.
#'   \item local - use the locfit package to fit a local regression
#'     of log dispersions over log base mean (normal scale means and dispersions
#'     are input and output for \code{\link{dispersionFunction}}). The points
#'     are weighted by normalized mean count in the local regression.
#'   \item mean - use the mean of gene-wise dispersion estimates.
#'   \item glmGamPoi - use the glmGamPoi package to fit the gene-wise dispersion, its trend
#'     and calculate the MAP based on the quasi-likelihood framework. The trend is 
#'     calculated using a local median regression.
#' }
#' @param maxit control parameter: maximum number of iterations to allow for convergence
#' @param useCR whether to use Cox-Reid correction - see McCarthy et al (2012)
#' @param weightThreshold threshold for subsetting the design matrix and GLM weights
#' for calculating the Cox-Reid correction
#' @param quiet whether to print messages at each step
#' @param modelMatrix an optional matrix which will be used for fitting the expected counts.
#' by default, the model matrix is constructed from \code{design(object)}
#' @param minmu lower bound on the estimated count for fitting gene-wise dispersion
#'
#' @return The DESeqDataSet passed as parameters, with the dispersion information
#' filled in as metadata columns, accessible via \code{mcols}, or the final dispersions
#' accessible via \code{\link{dispersions}}.
#'
#' @references \itemize{
#'   \item Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data.
#' Genome Biology 11 (2010) R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#'   \item McCarthy, DJ, Chen, Y, Smyth, GK: Differential expression analysis of multifactor RNA-Seq
#' experiments with respect to biological variation. Nucleic Acids Research 40 (2012), 4288-4297,
#' \url{http://dx.doi.org/10.1093/nar/gks042}
#'   \item Wu, H., Wang, C. & Wu, Z. A new shrinkage estimator for dispersion improves differential
#' expression detection in RNA-seq data. Biostatistics (2012).
#' \url{http://dx.doi.org/10.1093/biostatistics/kxs033}
#'   \item Ahlmann-Eltze, C., Huber, W. glmGamPoi: Fitting Gamma-Poisson Generalized Linear Models on Single Cell Count Data. Bioinformatics (2020).
#' \url{https://doi.org/10.1093/bioinformatics/btaa1009}
#' }
#'
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet()
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' head(dispersions(dds))
#'
#' @export
setMethod("estimateDispersions", signature(object="DESeqDataSet"),
          estimateDispersions.DESeqDataSet)


#' Show method for DESeqResults objects
#'
#' Prints out the information from the metadata columns
#' of the results object regarding the log2 fold changes
#' and p-values, then shows the DataFrame using the
#' standard method.
#' 
#' @rdname show
#'
#' @author Michael Love
#' 
#' @param object a DESeqResults object
#' 
#' @export
setMethod("show", signature(object="DESeqResults"), function(object) {
  cat(mcols(object)$description[ colnames(object) == "log2FoldChange"],"\n")
  cat(mcols(object)$description[ colnames(object) == "pvalue"],"\n")
  show(DataFrame(object))
})

#' Extract a matrix of model coefficients/standard errors
#'
#' \strong{Note:} results tables with log2 fold change, p-values, adjusted p-values, etc.
#' for each gene are best generated using the \code{\link{results}} function. The \code{coef}
#' function is designed for advanced users who wish to inspect all model coefficients at once.
#' 
#' Estimated model coefficients or estimated standard errors are provided in a matrix
#' form, number of genes by number of parameters, on the log2 scale.
#' The columns correspond to columns of the model matrix for final GLM fitting, i.e.,
#' \code{attr(dds, "modelMatrix")}.
#'
#' @param object a DESeqDataSet returned by \code{\link{DESeq}}, \code{\link{nbinomWaldTest}},
#' or \code{\link{nbinomLRT}}.
#' @param SE whether to give the standard errors instead of coefficients.
#' defaults to FALSE so that the coefficients are given.
#' @param ... additional arguments
#'
#' @rdname coef
#'
#' @author Michael Love
#' @importFrom stats coef
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' dds <- DESeq(dds)
#' coef(dds)[1,]
#' coef(dds, SE=TRUE)[1,]
#'
#' @method coef DESeqDataSet
#' @export
coef.DESeqDataSet  <- function(object, SE=FALSE, ...) {
  # Temporary hack for backward compatibility with "old" DESeqDataSet
  # objects. Remove once all serialized DESeqDataSet objects around have
  # been updated.
  if (!.hasSlot(object, "rowRanges"))
    object <- updateObject(object)
  resNms <- resultsNames(object)
  if (length(resNms) == 0) {
    stop("no coefficients have been generated yet, first call DESeq()")
  }
  if (!SE) {
    as.matrix(mcols(object,use.names=TRUE)[resNms])
  } else {
    as.matrix(mcols(object,use.names=TRUE)[paste0("SE_",resNms)])
  }
}

summary.DESeqResults <- function(object, alpha, ...) {
  sval <- "svalue" %in% names(object)
  if (sval) {
    test.col <- "svalue"
    test.col.name <- "s-value"
  } else {
    test.col <- "padj"
    test.col.name <- "adjusted p-value"
  }
  if (missing(alpha)) {
    if (sval) {
      alpha <- 0.005
    } else {
      if (is.null(metadata(object)$alpha)) {
        alpha <- 0.1
      } else {
        alpha <- metadata(object)$alpha
      }
    }
  }
  if (!is.null(metadata(object)$lfcThreshold)) {
    T <- metadata(object)$lfcThreshold
    pT <- sprintf("%.2f (up)    ", T)
    mT <- sprintf("%.2f (down) ", -T)
  } else {
    T <- 0
  }
  if (T == 0) {
    pT <- "0 (up)       "
    mT <- "0 (down)     "
  }
  cat("\n")
  notallzero <- sum(object$baseMean > 0)
  up <- sum(object[[test.col]] < alpha & object$log2FoldChange > T, na.rm=TRUE)
  down <- sum(object[[test.col]] < alpha & object$log2FoldChange < -T, na.rm=TRUE)
  if (!sval) {
    filt <- sum(!is.na(object$pvalue) & is.na(object$padj))
    outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
    if (is.null(metadata(object)$filterThreshold)) {
      ft <- 0
    } else {
      ft <- round(metadata(object)$filterThreshold)
    }
  }
  ihw <- !sval & "ihwResult" %in% names(metadata(object))
  
  printsig <- function(x) format(x, digits=2) 
  cat(paste("out of",notallzero,"with nonzero total read count\n"))
  cat(paste0(test.col.name," < ",alpha,"\n"))
  cat(paste0("LFC > ",pT,": ",up,", ",printsig(up/notallzero*100),"%\n"))
  cat(paste0("LFC < ",mT,": ",down,", ",printsig(down/notallzero*100),"%\n"))
  if (!sval) cat(paste0("outliers [1]       : ",outlier,", ",printsig(outlier/notallzero*100),"%\n"))
  if (!sval & !ihw) cat(paste0("low counts [2]     : ",filt,", ",printsig(filt/notallzero*100),"%\n"))
  if (!sval & !ihw) cat(paste0("(mean count < ",ft,")\n"))
  if (!sval) cat("[1] see 'cooksCutoff' argument of ?results\n")
  if (!sval & !ihw) cat("[2] see 'independentFiltering' argument of ?results\n")
  if (ihw) cat("see metadata(res)$ihwResult on hypothesis weighting\n")
  cat("\n")
}

#' Summarize DESeq results
#'
#' Print a summary of the results from a DESeq analysis.
#'
#' @param object a \code{\link{DESeqResults}} object
#' @param alpha the adjusted p-value cutoff. If not set, this
#' defaults to the \code{alpha} argument which was used in
#' \code{\link{results}} to set the target FDR for independent
#' filtering, or if independent filtering was not performed,
#' to 0.1.
#' @param ... additional arguments
#'
#' @rdname summary
#'
#' @author Michael Love
#'  
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' summary(res)
#'
#' @method summary DESeqResults
#' @export
setMethod("summary", signature(object="DESeqResults"), summary.DESeqResults)

#' Accessors for the 'priorInfo' slot of a DESeqResults object.
#' 
#' The priorInfo slot contains details about the prior on log fold changes
#' 
#' @rdname priorInfo
#' 
#' @param object a \code{DESeqResults} object
#' @param value a \code{list}
#' @param ... additional arguments
#'
#' @export
setMethod("priorInfo", signature(object="DESeqResults"),
          function(object) object@priorInfo)

#' @rdname priorInfo
#' @export
setReplaceMethod("priorInfo", signature(object="DESeqResults", value="list"),
                 function(object, value) {
                   object@priorInfo <- value
                   object
                 })
