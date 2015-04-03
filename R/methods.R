counts.DESeqDataSet <- function(object, normalized=FALSE, replaced=FALSE) {
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
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,DESeqDataSet-method counts<-,DESeqDataSet,matrix-method
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

#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
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
#' @docType methods
#' @name design
#' @rdname design
#' @aliases design design,DESeqDataSet-method design<-,DESeqDataSet,formula-method
#' @param object a \code{DESeqDataSet} object
#' @param value a \code{formula} used for estimating dispersion
#' and fitting Negative Binomial GLMs
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' design(dds) <- formula(~ 1)
#'
#' @export
setMethod("design", signature(object="DESeqDataSet"), design.DESeqDataSet)

#' @name design
#' @rdname design
#' @exportMethod "design<-"
setReplaceMethod("design", signature(object="DESeqDataSet", value="formula"),
                 function( object, value ) {  
                   object@design <- value
                   validObject(object)
                   object
                 })

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
#' @docType methods
#' @name dispersionFunction
#' @rdname dispersionFunction
#' @aliases dispersionFunction dispersionFunction,DESeqDataSet-method dispersionFunction<-,DESeqDataSet,function-method
#' @param object a \code{DESeqDataSet} object.
#' @param value a \code{function}
#' @param estimateVar whether to estimate the variance of dispersion residuals.
#' setting to FALSE is needed, e.g. within \code{estimateDispersionsMAP} when
#' called on a subset of the full dataset in parallel execution.
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

#' @name dispersionFunction
#' @rdname dispersionFunction
#' @exportMethod "dispersionFunction<-"
setReplaceMethod("dispersionFunction",
                 signature(object="DESeqDataSet", value="function"),
                 function(object, value, estimateVar=TRUE) {

                   if (estimateVar) {
                     if (is.null(mcols(object)$baseMean) | is.null(mcols(object)$allZero)) {
                       object <- getBaseMeansAndVariances(object)
                     }

                     if (!is.null(mcols(object)$dispFit)) {
                       message("found already estimated fitted dispersions, removing these")
                       mcols(object) <- mcols(object)[,!names(mcols(object)) == "dispFit",drop=FALSE]
                     }
                     
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
                   }
                   
                   object@dispersionFunction <- value   
                   validObject(object)
                   object
                 })

dispersions.DESeqDataSet <- function(object) mcols(object)$dispersion

#' Accessor functions for the dispersion estimates in a DESeqDataSet
#' object.
#' 
#' The dispersions for each row of the DESeqDataSet.  Generally,
#' these are set by \code{\link{estimateDispersions}}.
#' 
#' @docType methods
#' @name dispersions
#' @rdname dispersions
#' @aliases dispersions dispersions,DESeqDataSet-method dispersions<-,DESeqDataSet,numeric-method
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

#' @name dispersions
#' @rdname dispersions
#' @exportMethod "dispersions<-"
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
#' See \code{\link{DESeq}} for a description of the use of size factors. If gene-specific normalization
#' is desired for each sample, use \code{\link{normalizationFactors}}.
#' 
#' @docType methods
#' @name sizeFactors
#' @rdname sizeFactors
#' @aliases sizeFactors sizeFactors,DESeqDataSet-method sizeFactors<-,DESeqDataSet,numeric-method
#' @param object a \code{DESeqDataSet} object.
#' @param value a numeric vector, one size factor for each column in the count
#' data.
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactors}}
#'
#' @export
setMethod("sizeFactors", signature(object="DESeqDataSet"),
          sizeFactors.DESeqDataSet)

#' @name sizeFactors
#' @rdname sizeFactors
#' @exportMethod "sizeFactors<-"
setReplaceMethod("sizeFactors", signature(object="DESeqDataSet", value="numeric"),
                 function( object, value ) {
                   if (any(value <= 0)) {
                     stop("size factors must be positive")
                   }
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

#' @name normalizationFactors
#' @rdname normalizationFactors
#' @exportMethod "normalizationFactors<-"
setReplaceMethod("normalizationFactors", signature(object="DESeqDataSet", value="matrix"),
                 function(object, value) {
                   if (any(value <= 0)) {
                     stop("normalization factors must be positive")
                   }
                   # enforce same dimnames
                   dimnames(value) <- dimnames(object)
                   assays(object)[["normalizationFactors"]] <- value
                   validObject( object )
                   object
                 })

estimateSizeFactors.DESeqDataSet <- function(object, type=c("ratio","iterate"),
                                             locfunc=stats::median, geoMeans, controlGenes, normMatrix) {
  type <- match.arg(type, c("ratio","iterate"))
  object <- sanitizeColData(object)
  if (type == "iterate") {
    sizeFactors(object) <- estimateSizeFactorsIterate(object)
  } else {
    if (missing(normMatrix)) {
      sizeFactors(object) <- estimateSizeFactorsForMatrix(counts(object), locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
    } else {
      normalizationFactors(object) <- estimateNormFactors(counts(object), normMatrix=normMatrix,
                                                          locfunc=locfunc,
                                                          geoMeans=geoMeans,
                                                          controlGenes=controlGenes)
      message("adding normalization factors which account for library size")
    }
  }
  object
}

#' Estimate the size factors for a DESeqDataSet
#' 
#' This function estimates the size factors using the
#' "median ratio method" described by Equation 5 in Anders and Huber (2010).
#' The estimated size factors can be accessed using \code{\link{sizeFactors}}.
#' Alternative library size estimators can also be supplied
#' using \code{\link{sizeFactors}}.
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
#' @docType methods
#' @name estimateSizeFactors
#' @rdname estimateSizeFactors
#' @aliases estimateSizeFactors estimateSizeFactors,DESeqDataSet-method
#' 
#' @param object a DESeqDataSet
#' @param type either "ratio" or "iterate". "ratio" uses the standard
#' median ratio method introduced in DESeq. The size factor is the
#' median ratio of the sample over a pseudosample: for each gene, the geometric mean
#' of all samples. "iterate" offers an alternative estimator, which can be
#' used even when all genes contain a sample with a zero. This estimator
#' iterates between estimating the dispersion with a design of ~1, and
#' finding a size factor vector by numerically optimizing the likelihood
#' of the ~1 model.
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used. However, especially for low counts, the
#' \code{\link[genefilter]{shorth}} function from the genefilter package may give better results.
#' @param geoMeans by default this is not provided and the
#' geometric means of the counts are calculated within the function.
#' A vector of geometric means from another count matrix can be provided
#' for a "frozen" size factor calculation
#' @param controlGenes optional, numeric or logical index vector specifying those genes to
#' use for size factor estimation (e.g. housekeeping or spike-in genes)
#' @param normMatrix optional, a matrix of normalization factors which do not
#' control for library size (e.g. average transcript length of genes for each
#' sample). Providing \code{normMatrix} will estimate size factors on the
#' count matrix divided by \code{normMatrix} and store the product of the
#' size factors and \code{normMatrix} as \code{\link{normalizationFactors}}.
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
#' Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 2010, 11:106. \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
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

estimateDispersions.DESeqDataSet <- function(object, fitType=c("parametric","local","mean"),
                                             maxit=100, quiet=FALSE, modelMatrix=NULL) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
  }
  # size factors could have slipped in to colData from a previous run
  if (!is.null(sizeFactors(object))) {
    if (!is.numeric(sizeFactors(object))) {
      stop("the sizeFactor column in colData is not numeric.
these column could have come in during colData import")
    }
    if (any(is.na(sizeFactors(object)))) {
      stop("the sizeFactor column in colData contains NA.
these column could have come in during colData import")
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
  fitType <- match.arg(fitType, choices=c("parametric","local","mean"))
  
  noReps <- checkForExperimentalReplicates(object, modelMatrix)
  if (noReps) {
    designIn <- design(object)
    design(object) <- formula(~ 1)
  }
  
  if (!quiet) message("gene-wise dispersion estimates")
  object <- estimateDispersionsGeneEst(object, maxit=maxit, quiet=quiet, modelMatrix=modelMatrix)
  if (!quiet) message("mean-dispersion relationship")
  object <- estimateDispersionsFit(object, fitType=fitType, quiet=quiet)
  if (!quiet) message("final dispersion estimates")
  object <- estimateDispersionsMAP(object, maxit=maxit, quiet=quiet, modelMatrix=modelMatrix)

  # replace the previous design
  if (noReps) design(object) <- designIn
  
  return(object)
}

checkForExperimentalReplicates <- function(object, modelMatrix) {
  noReps <- if (is.null(modelMatrix)) {
    mmtest <- model.matrix(design(object), data=as.data.frame(colData(object)))
    nrow(mmtest) == ncol(mmtest)
  } else {
    nrow(modelMatrix) == ncol(modelMatrix)
  }
  if (noReps) {
    if (!is.null(modelMatrix)) stop("same number of samples and coefficients to fit with supplied model matrix")
    warning("same number of samples and coefficients to fit,
  estimating dispersion by treating samples as replicates.
  read the ?DESeq section on 'Experiments without replicates'")
  }
  noReps
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
#' \code{estimateDispersions} checks for the case of an analysis
#' with as many samples as the number of coefficients to fit,
#' and will temporarily substitute a design formula \code{~ 1} for the
#' purposes of dispersion estimation.  This treats the samples as 
#' replicates for the purpose of dispersion estimation. As mentioned in the DESeq paper:
#' "While one may not want to draw strong conclusions from such an analysis,
#' it may still be useful for exploration and hypothesis generation."
#'
#' The lower-level functions called by \code{estimateDispersions} are:
#' \code{\link{estimateDispersionsGeneEst}},
#' \code{\link{estimateDispersionsFit}}, and
#' \code{\link{estimateDispersionsMAP}}.
#' 
#' @docType methods
#' @name estimateDispersions
#' @rdname estimateDispersions
#' @aliases estimateDispersions estimateDispersions,DESeqDataSet-method
#' @param object a DESeqDataSet
#' @param fitType either "parametric", "local", or "mean"
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
#' }
#' @param maxit control parameter: maximum number of iterations to allow for convergence
#' @param quiet whether to print messages at each step
#' @param modelMatrix an optional matrix which will be used for fitting the expected counts.
#' by default, the model matrix is constructed from \code{design(object)}
#'
#' @return The DESeqDataSet passed as parameters, with the dispersion information
#' filled in as metadata columns, accessible via \code{mcols}, or the final dispersions
#' accessible via \code{\link{dispersions}}.
#'
#' @references \itemize{
#'   \item Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 11 (2010) R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#'   \item McCarthy, DJ, Chen, Y, Smyth, GK: Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. Nucleic Acids Research 40 (2012), 4288-4297, \url{http://dx.doi.org/10.1093/nar/gks042}
#'   \item Wu, H., Wang, C. & Wu, Z. A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data. Biostatistics (2012). \url{http://dx.doi.org/10.1093/biostatistics/kxs033}
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
#' @docType methods
#' @name show
#' @rdname show
#' @aliases show show,DESeqResults-method
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
#' @docType methods
#' @name coef
#' @rdname coef
#' @aliases coef coef.DESeqDataSet
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
#' @export
coef.DESeqDataSet  <- function(object, SE=FALSE, ...) {
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

#' Summarize DESeq results
#'
#' Print a summary of the results from a DESeq analysis.
#'
#' @usage
#' \method{summary}{DESeqResults}(object, alpha=.1, \dots)
#' 
#' @param object a \code{\link{DESeqResults}} object
#' @param alpha the adjusted p-value cutoff
#' @param ... additional arguments
#'
#' @docType methods
#' @name summary
#' @rdname summary
#' @aliases summary summary.DESeqResults
#' @author Michael Love
#'  
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(m=4)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' summary(res)
#'
#' @export
summary.DESeqResults <- function(object, alpha=.1, ...) {
  cat("\n")
  notallzero <- sum(object$baseMean > 0)
  up <- sum(object$padj < alpha & object$log2FoldChange > 0, na.rm=TRUE)
  down <- sum(object$padj < alpha & object$log2FoldChange < 0, na.rm=TRUE)
  filt <- sum(!is.na(object$pvalue) & is.na(object$padj))
  outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
  ft <- if (is.null(attr(object, "filterThreshold"))) {
    0
  } else {
    round(attr(object,"filterThreshold"), 1)
  }
  printsig <- function(x) format(x, digits=2) 
  cat("out of",notallzero,"with nonzero total read count\n")
  cat(paste0("adjusted p-value < ",alpha,"\n"))
  cat(paste0("LFC > 0 (up)     : ",up,", ",printsig(up/notallzero*100),"% \n"))
  cat(paste0("LFC < 0 (down)   : ",down,", ",printsig(down/notallzero*100),"% \n"))
  cat(paste0("outliers [1]     : ",outlier,", ",printsig(outlier/notallzero*100),"% \n"))
  cat(paste0("low counts [2]   : ",filt,", ",printsig(filt/notallzero*100),"% \n"))
  cat(paste0("(mean count < ",ft,")\n"))
  cat("[1] see 'cooksCutoff' argument of ?results\n")
  cat("[2] see 'independentFiltering' argument of ?results\n")
  cat("\n")
}
