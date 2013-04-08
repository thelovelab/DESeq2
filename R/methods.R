#' Accessors for the 'counts' slot of a DESeqSummarizedExperiment object.
#' 
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample. 
#'
#' @usage
#' \S4method{counts}{DESeqSummarizedExperiment}(object,normalized=FALSE)
#'
#' \S4method{counts}{DESeqSummarizedExperiment,matrix}(object)<-value
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,DESeqSummarizedExperiment-method counts<-,DESeqSummarizedExperiment,matrix-method
#' @param object a \code{DESeqSummarizedExperiment} object.
#' @param normalized logical indicating whether or not to divide the counts by
#' the size factors or normalization factors before returning
#' (normalization factors always preempt size factors)
#' @param value an integer matrix
#' @author Simon Anders
#' @seealso \code{\link{sizeFactors}}, \code{\link{normalizationFactors}}
#'
#' @examples
#' 
#' dse <- makeExampleDESeqSummarizedExperiment()
#' head(counts(dse))
#'
counts.DESeqSummarizedExperiment <- function(object, normalized=FALSE) {
            if(!normalized) {
              return(assays(object)[["counts"]])
            } else {
              if (!is.null(normalizationFactors(object))) {
                return( assays(object)[["counts"]]/normalizationFactors(object) )
              } else if (is.null(sizeFactors(object)) | any(is.na(sizeFactors(object)))) {
                stop("first calculate size factors, add normalizationFactors, or set normalized=FALSE")
              } else {
                return( t( t( assays(object)[["counts"]] ) / sizeFactors(object) ) )
              }
            }
          }

#' @rdname counts
#' @export
setMethod("counts", signature(object="DESeqSummarizedExperiment"), counts.DESeqSummarizedExperiment)

#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object="DESeqSummarizedExperiment", value="matrix"),
  function( object, value ) {
   assays(object)[["counts"]] <- value
   validObject(object)
   object
})   


#' Accessors for the 'design' slot of a DESeqSummarizedExperiment object.
#' 
#' Accessors for the 'design' slot of a DESeqSummarizedExperiment object.
#' 
#' @usage
#' \S4method{design}{DESeqSummarizedExperiment}(object)
#'
#' \S4method{design}{DESeqSummarizedExperiment,formula}(object)<-value
#'
#' @docType methods
#' @name design
#' @rdname design
#' @aliases design design,DESeqSummarizedExperiment-method design<-,DESeqSummarizedExperiment,formula-method
#' @param object a \code{DESeqSummarizedExperiment} object
#' @param value a \code{formula} used for estimating dispersion
#' and fitting negative binomial GLMs
#' @examples
#'
#' dse <- makeExampleDESeqSummarizedExperiment()
#' design(dse) <- formula(~ 1)
#'
design.DESeqSummarizedExperiment <- function(object) object@design

#' @rdname design
#' @export
setMethod("design", signature(object="DESeqSummarizedExperiment"), design.DESeqSummarizedExperiment)

#' @name design
#' @rdname design
#' @exportMethod "design<-"
setReplaceMethod("design", signature(object="DESeqSummarizedExperiment", value="formula"),
  function( object, value ) {  
    object@design <- value
    validObject(object)
    object
})


#' Accessors for the 'dispersionFunction' slot of a DESeqSummarizedExperiment object.
#'
#' The dispersion function is calculated by \code{\link{estimateDispersions}} and
#' used by \code{\link{varianceStabilizingTransformation}}.  Parametric dispersion
#' fits store the coefficients of the fit as attributes in this slot.
#'
#' @usage
#' \S4method{dispersionFunction}{DESeqSummarizedExperiment}(object)
#'
#' \S4method{dispersionFunction}{DESeqSummarizedExperiment,function}(object)<-value
#'
#' @docType methods
#' @name dispersionFunction
#' @rdname dispersionFunction
#' @aliases dispersionFunction dispersionFunction,DESeqSummarizedExperiment-method dispersionFunction<-,DESeqSummarizedExperiment,function-method
#' @param object a \code{DESeqSummarizedExperiment} object.
#' @param value a \code{function}
#' @examples
#'
#' example("estimateDispersions")
#' dispersionFunction(dse)
#'
dispersionFunction.DESeqSummarizedExperiment <- function(object) object@dispersionFunction

#' @rdname dispersionFunction
#' @export
setMethod("dispersionFunction", signature(object="DESeqSummarizedExperiment"),
          dispersionFunction.DESeqSummarizedExperiment)
setReplaceMethod("dispersionFunction", signature(object="DESeqSummarizedExperiment", value="function"),
  function( object, value ) {  
    object@dispersionFunction <- value
    validObject(object)
    object
})   

#' Accessor functions for the dispersion estimates in a DESeqSummarizedExperiment
#' object.
#' 
#' The dispersions for each row of the DESeqSummarizedExperiment.  Generally,
#' these should be set only by \code{\link{estimateDispersions}}.
#' 
#' @usage
#' \S4method{dispersions}{DESeqSummarizedExperiment}(object)
#'
#' \S4method{dispersions}{DESeqSummarizedExperiment,numeric}(object)<-value
#'
#' @docType methods
#' @name dispersions
#' @rdname dispersions
#' @aliases dispersions dispersions,DESeqSummarizedExperiment-method dispersions<-,DESeqSummarizedExperiment,numeric-method
#' @param object a \code{DESeqSummarizedExperiment} object.
#' @param value the dispersions to use for the negative binomial modeling
#'
#' @author Simon Anders
#' @seealso \code{\link{estimateDispersions}}
#' @examples
#'
#' example("estimateDispersions")
#' dispersions(dse)
#'
dispersions.DESeqSummarizedExperiment <- function(object) mcols(object)$dispersion

#' @rdname dispersions
#' @export
setMethod("dispersions", signature(object="DESeqSummarizedExperiment"),
          dispersions.DESeqSummarizedExperiment)
setReplaceMethod("dispersions", signature(object="DESeqSummarizedExperiment", value="numeric"),
                 function(object, value) {
                   mcols(object)$dispersion <- value
                   validObject( object )
                   object
                 })



#' Accessor functions for the 'sizeFactors' information in a DESeqSummarizedExperiment
#' object.
#' 
#' The sizeFactors vector assigns to each column of the count matrix a value, the
#' size factor, such that count values in the columns can be brought to a common
#' scale by dividing by the corresponding size factor.  See \code{\link{DESeq}}
#' for a description of the use of size factors. If gene-specific normalization
#' is desired for each sample, use \code{\link{normalizationFactors}}.
#' 
#' @usage
#' \S4method{sizeFactors}{DESeqSummarizedExperiment}(object)
#'
#' \S4method{sizeFactors}{DESeqSummarizedExperiment,numeric}(object)<-value
#'
#' @docType methods
#' @name sizeFactors
#' @rdname sizeFactors
#' @aliases sizeFactors sizeFactors,DESeqSummarizedExperiment-method sizeFactors<-,DESeqSummarizedExperiment,numeric-method
#' @param object a \code{DESeqSummarizedExperiment} object.
#' @param value a numeric vector, one size factor for each column in the count
#' data.
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactors}}
#' @examples
#' 
#' dse <- makeExampleDESeqSummarizedExperiment()
#' dse <- estimateSizeFactors( dse )
#' sizeFactors(dse)
sizeFactors.DESeqSummarizedExperiment <- function(object) {
  if (!"sizeFactor" %in% names(colData(object))) return(NULL)
  sf <- colData(object)$sizeFactor
  names( sf ) <- colnames( object )
  sf
}

#' @rdname sizeFactors
#' @export
setMethod("sizeFactors", signature(object="DESeqSummarizedExperiment"),
          sizeFactors.DESeqSummarizedExperiment)

#' @name sizeFactors
#' @rdname sizeFactors
#' @exportMethod "sizeFactors<-"
setReplaceMethod("sizeFactors", signature(object="DESeqSummarizedExperiment", value="numeric"),
                 function( object, value ) {
                   if (any(value <= 0)) {
                     stop("size factors must be positive")
                   }
                   # have to make sure to remove sizeFactor which might be
                   # coming from a previous CountDataSet
                   colData(object)$sizeFactor <- value
                   idx <- which(colnames(colData(object)) == "sizeFactor")
                   metaDataFrame <- DataFrame(type="intermediate",
                                              description="a scaling factor for columns")
                   mcols(colData(object))[idx,] <- metaDataFrame
                   validObject( object )
                   object
                 }) 

#' Accessor functions for the normalization factors in a DESeqSummarizedExperiment
#' object.
#'
#' Gene-specific normalization factors for each sample can be provided as a matrix,
#' which will preempt \code{\link{sizeFactors}}. In some experiments, counts for each
#' sample have varying dependence on covariates, e.g. on GC-content for sequencing
#' data run on different days, and in this case it makes sense to provide
#' gene-specific factors for each sample rather than a single size factor.
#'
#' Normalization factors alter the model of \code{\link{DESeq}} in the following way, for
#' counts \eqn{K_{ij}}{K_ij} for gene i and sample j:
#'
#' \deqn{K_{ij} \sim \textrm{NB}(NF_{ij} \mu_{ij}, \alpha_i)}{K_ij ~ NB(NF_ij * mu_ij, alpha_i)}
#'
#' @note Normalization factors are on the scale of the counts (similar to \code{\link{sizeFactors}})
#' and unlike offsets, which are typically on the scale of the predictors (in this case, log counts).
#' Normalization factors should include size factor normalization and should have
#' a mean around 1, as is the case with size factors.
#'
#' @usage
#' \S4method{normalizationFactors}{DESeqSummarizedExperiment}(object)
#'
#' \S4method{normalizationFactors}{DESeqSummarizedExperiment,matrix}(object)<-value
#'
#' @docType methods
#' @name normalizationFactors
#' @rdname normalizationFactors
#' @aliases normalizationFactors normalizationFactors,DESeqSummarizedExperiment-method normalizationFactors<-,DESeqSummarizedExperiment,matrix-method
#' @param object a \code{DESeqSummarizedExperiment} object.
#' @param value the matrix of normalization factors
#' @examples
#'
#' dse <- makeExampleDESeqSummarizedExperiment()
#' normFactors <- matrix(runif(nrow(dse)*ncol(dse),0.5,1.5),
#'                       ncol=ncol(dse),nrow=nrow(dse))
#' normalizationFactors(dse) <- normFactors
#' dse <- estimateDispersions(dse)
#' dse <- nbinomWaldTest(dse)
#'
normalizationFactors.DESeqSummarizedExperiment <- function(object) {
  if (!"normalizationFactors" %in% names(assays(object))) return(NULL)
  assays(object)[["normalizationFactors"]]
}

#' @rdname normalizationFactors
#' @export
setMethod("normalizationFactors", signature(object="DESeqSummarizedExperiment"),
          normalizationFactors.DESeqSummarizedExperiment)
setReplaceMethod("normalizationFactors", signature(object="DESeqSummarizedExperiment", value="matrix"),
                 function(object, value) {
                   if (any(value <= 0)) {
                     stop("normalization factors must be positive")
                   }
                   if (mean(value) > 10) {
                     warning("the replacement matrix has a mean above 10, normalization factors should have an overall mean near 1 as is the case with size factors")
                   }
                   assays(object)[["normalizationFactors"]] <- value
                   validObject( object )
                   object
                 })


#' Estimate the size factors for a DESeqSummarizedExperiment
#' 
#' Estimate the size factors for a DESeqSummarizedExperiment
#' 
#' This function estimates the size factors and stores the information
#' which can be accessed using \code{\link{sizeFactors}}
#'
#' Typically, the function is called with the idiom:
#'
#' \code{dse <- estimateSizeFactors(dse)}
#'
#' See \code{\link{DESeq}} for a description of the use of size factors in the GLM.
#' You need to call this function after \code{\link{DESeqSummarizedExperiment}}
#' unless you have manually specified \code{\link{sizeFactors}}.
#' Alternatively, gene-specific normalization factors for each sample can be provided using
#' \code{\link{normalizationFactors}} which will always preempt \code{\link{sizeFactors}}
#' in calculations.
#'
#' Internally, the function calls \code{\link{estimateSizeFactorsForMatrix}}, 
#' which provides more details on the calculation.
#'
#' @usage
#' \S4method{estimateSizeFactors}{DESeqSummarizedExperiment}(object,locfunc=median)
#'
#' @docType methods
#' @name estimateSizeFactors
#' @rdname estimateSizeFactors
#' @aliases estimateSizeFactors estimateSizeFactors,DESeqSummarizedExperiment-method
#' @param object a DESeqSummarizedExperiment
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used. However, especially for low counts, the
#' \code{\link[genefilter]{shorth}} function from the genefilter package may give better results.
#' @return The DESeqSummarizedExperiment passed as parameters, with the size factors filled
#' in.
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactorsForMatrix}}
#' 
#' @examples
#' 
#' dse <- makeExampleDESeqSummarizedExperiment()
#' dse <- estimateSizeFactors( dse )
#' sizeFactors( dse )
estimateSizeFactors.DESeqSummarizedExperiment <- function(object, locfunc=median) {
  sizeFactors(object) <- estimateSizeFactorsForMatrix( counts(object), locfunc )
  object
}
  
#' @rdname estimateSizeFactors
#' @export
setMethod("estimateSizeFactors", signature(object="DESeqSummarizedExperiment"),
          estimateSizeFactors.DESeqSummarizedExperiment)


#' Estimate the dispersions for a DESeqSummarizedExperiment
#' 
#' This function obtains dispersion estimates for negative binomial distributed data.
#'
#' Typically the function is called with the idiom:
#'
#' \code{dse <- estimateDispersions(dse)}
#'
#' The fitting proceeds as follows: for each gene, an estimate of the dispersion
#' is found which maximizes the Cox Reid-adjusted profile likelihood
#' (the methods of Cox Reid-adjusted profile likelihood maximization for
#' estimation of dispersion in RNA-Seq data were developed by McCarthy,
#' et al. (2012), first implemented in the edgeR package in 2010);
#' a dispersion-mean relationship is fit to the maximum likelihood estimates;
#' a normal prior is determined for the log dispersion estimates centered
#' on the predicted value from the fit
#' with variance equal to the difference between the observed variance of the
#' log dispersion estimates and the expected sampling variance;
#' finally maximum a posteriori dispersion estimates are returned.
#' This final dispersion parameter is used in subsequent tests.
#' The final dispersion estimates can be accessed from an object using \code{\link{dispersions}}.
#' The fitted dispersion-mean relationship is also used in
#' \code{\link{varianceStabilizingTransformation}}.
#'
#' The log normal prior on the dispersion parameter has been proposed
#' by Wu, et al. (2012) and is also implemented in the DSS package.
#'
#' \code{estimateDispersions} checks for the case of a 1-vs-1
#' comparison, and will temporarily substitute a design formula \code{~ 1} for the
#' purposes of dispersion estimation.  This treats the two samples as 
#' replicates for the purpose of dispersion estimation. As mentioned in the DESeq paper:
#' "While one may not want to draw strong conclusions from such an analysis,
#' it may still be useful for exploration and hypothesis generation."
#'
#' The lower-level functions called by \code{estimateDispersions} are:
#' \code{\link{estimateDispersionsGeneEst}},
#' \code{\link{estimateDispersionsFit}}, and
#' \code{\link{estimateDispersionsMAP}}.
#' 
#' @usage
#' \S4method{estimateDispersions}{DESeqSummarizedExperiment}(object,fitType=c("parametric","local","geoMean"))
#'
#' @docType methods
#' @name estimateDispersions
#' @rdname estimateDispersions
#' @aliases estimateDispersions estimateDispersions,DESeqSummarizedExperiment-method
#' @param object a DESeqSummarizedExperiment
#' @param fitType either "parametric", "local", or "geoMean"
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
#'   \item geoMean - use the geometric mean of gene-wise dispersion estimates.
#'     This is motivated by an assumption that the log normal prior contributes
#'     more to the observed spread of dispersion estimates than the sampling
#'     variance.
#' }
#'
#' @return The DESeqSummarizedExperiment passed as parameters, with the dispersion information
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
#' dse <- makeExampleDESeqSummarizedExperiment()
#' dse <- estimateSizeFactors(dse)
#' dse <- estimateDispersions(dse)
#' head(dispersions(dse))
#'
estimateDispersions.DESeqSummarizedExperiment <- function(object, fitType=c("parametric","local","geoMean")) {
  if (is.null(sizeFactors(object)) & is.null(normalizationFactors(object))) {
    stop("first call estimateSizeFactors or provide a normalizationFactor matrix before estimateDispersions")
  }
  if (!is.null(dispersions(object))) {
    message("you had estimated dispersions, replacing these")
    mcols(object) <- mcols(object)[,!(mcols(mcols(object))$type %in% c("intermediate","results"))]
  }
  fitType <- match.arg(fitType)
  
  # if trying to call differential expression for a model
  # with as many samples as columns which are factors, e.g.
  # 2 samples and 2 groups,
  # we supply a design formula of ~ 1 for dispersion estimation
  modelMatrix <- model.matrix(design(object), data=as.data.frame(colData(object)))  
  noReps <- nrow(modelMatrix) == sum(apply(modelMatrix, 2, function(x) all(x %in% 0:1)))
  if (noReps) {
    message("same number of samples and factor variables, estimating dispersion by treating samples as replicates")
    designIn <- design(object)
    design(object) <- formula(~ 1)
  }
  
  message("gene-wise dispersion estimates")
  object <- estimateDispersionsGeneEst(object)
  message("mean-dispersion relationship")
  object <- estimateDispersionsFit(object, fitType=fitType)
  message("final dispersion estimates")
  object <- estimateDispersionsMAP(object)

  # replace the previous design
  if (noReps) design(object) <- designIn
  
  return(object)
}

#' @rdname estimateDispersions
#' @export
setMethod("estimateDispersions", signature(object="DESeqSummarizedExperiment"),
          estimateDispersions.DESeqSummarizedExperiment)

