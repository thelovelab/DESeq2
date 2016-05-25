############################################################
#
# DESeq2 organization of R files
#
# core.R ......... most of the statistical code (example call below)
# methods.R ...... the S4 methods (estimateSizeFactors, etc.)
# AllClasses.R ... class definitions and object constructors
# AllGenerics.R .. the generics defined in DESeq2
# results.R ...... results() function and helpers
# plots.R ........ all plotting functions
# helper.R ....... collapseReplicates, fpkm, fpm, DESeqParallel
# expanded.R ..... helpers for dealing with expanded model matrices
# RcppExports.R .. the R wrappers for the C++ functions
# rlogTransformation.R
# varianceStabilizingTransformation.R
#
# general outline of the internal function calls.
# note: not all of these functions are exported.
#
# DESeq
# |- estimateSizeFactors
#    |- estimateSizeFactorsForMatrix
# |- estimateDispersions
#    |- estimateDispersionsGeneEst
#       |- fitNbinomGLMs
#          |- fitBeta (C++)
#       |- fitDisp (C++)
#    |- estimateDispersionsFit
#    |- estimateDispersionsMAP
#       |- estimateDispersionPriorVar
#       |- fitDisp (C++)
# |- nbinomWaldTest
#    |- fitGLMsWithPrior
#       |- fitNbinomGLMs
#          |- fitBeta (C++)
#       |- estimateBetaPriorVar
#       |- fitNbinomGLMs
#          |- fitBeta (C++)
#
############################################################

#' DESeq2 package for differential analysis of count data
#' 
#' The main functions for differential analysis are \code{\link{DESeq}} and
#' \code{\link{results}}. See the examples at \code{\link{DESeq}} for basic analysis steps.
#' Two transformations offered for count data are
#' the "regularized logarithm", \code{\link{rlog}},
#' and \code{\link{varianceStabilizingTransformation}}.
#' For more detailed information on usage, see the package vignette, by typing
#' \code{vignette("DESeq2")}, or the workflow linked to on the first page
#' of the vignette. All support questions should be posted to the Bioconductor
#' support site: \url{http://support.bioconductor.org}.
#'
#' @references
#'
#' DESeq2 reference:
#' 
#' Michael I Love, Wolfgang Huber, Simon Anders: Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 2014, 15:550. \url{http://dx.doi.org/10.1186/s13059-014-0550-8}
#'
#' DESeq reference:
#' 
#' Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 2010, 11:106. \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#'
#' @author Michael Love, Wolfgang Huber, Simon Anders
#' 
#' @docType package
#' @name DESeq2-package
#' @aliases DESeq2-package
#' @keywords package
NULL

#' Differential expression analysis based on the Negative Binomial (a.k.a. Gamma-Poisson) distribution
#'
#' This function performs a default analysis through the steps:
#' \enumerate{
#' \item estimation of size factors: \code{\link{estimateSizeFactors}}
#' \item estimation of dispersion: \code{\link{estimateDispersions}}
#' \item Negative Binomial GLM fitting and Wald statistics: \code{\link{nbinomWaldTest}}
#' }
#' For complete details on each step, see the manual pages of the respective
#' functions. After the \code{DESeq} function returns a DESeqDataSet object,
#' results tables (log2 fold changes and p-values) can be generated
#' using the \code{\link{results}} function. See the manual page
#' for \code{\link{results}} for information on independent filtering and
#' p-value adjustment for multiple test correction.
#'
#' The differential expression analysis uses a generalized linear model of the form:
#'
#' \deqn{ K_{ij} \sim \textrm{NB}( \mu_{ij}, \alpha_i) }{ K_ij ~ NB(mu_ij, alpha_i) }
#' \deqn{ \mu_{ij} = s_j q_{ij} }{ mu_ij = s_j q_ij }
#' \deqn{ \log_2(q_{ij}) = x_{j.} \beta_i }{ log2(q_ij) = x_j. beta_i }
#'
#' where counts \eqn{K_{ij}}{K_ij} for gene i, sample j are modeled using
#' a Negative Binomial distribution with fitted mean \eqn{\mu_{ij}}{mu_ij}
#' and a gene-specific dispersion parameter \eqn{\alpha_i}{alpha_i}.
#' The fitted mean is composed of a sample-specific size factor
#' \eqn{s_j}{s_j} and a parameter \eqn{q_{ij}}{q_ij} proportional to the
#' expected true concentration of fragments for sample j.
#' The coefficients \eqn{\beta_i}{beta_i} give the log2 fold changes for gene i for each
#' column of the model matrix \eqn{X}{X}.
#' The sample-specific size factors can be replaced by
#' gene-specific normalization factors for each sample using
#' \code{\link{normalizationFactors}}.
#'
#' For details on the fitting of the log2 fold changes and calculation of p-values,
#' see \code{\link{nbinomWaldTest}} if using \code{test="Wald"},
#' or \code{\link{nbinomLRT}} if using \code{test="LRT"}.
#'
#' Experiments without replicates do not allow for estimation of the dispersion
#' of counts around the expected value for each group, which is critical for
#' differential expression analysis. If an experimental design is
#' supplied which does not contain the necessary degrees of freedom for differential
#' analysis, \code{DESeq} will provide a message to the user and follow
#' the strategy outlined in Anders and Huber (2010)
#' under the section 'Working without replicates', wherein all the samples
#' are considered as replicates of a single group for the estimation of dispersion.
#' As noted in the reference above: "Some overestimation of the variance
#' may be expected, which will make that approach conservative."
#' Furthermore, "while one may not want to draw strong conclusions from such an analysis,
#' it may still be useful for exploration and hypothesis generation."
#'
#' The argument \code{minReplicatesForReplace} is used to decide which samples
#' are eligible for automatic replacement in the case of extreme Cook's distance.
#' By default, \code{DESeq} will replace outliers if the Cook's distance is
#' large for a sample which has 7 or more replicates (including itself).
#' This replacement is performed by the \code{\link{replaceOutliers}}
#' function. This default behavior helps to prevent filtering genes
#' based on Cook's distance when there are many degrees of freedom.
#' See \code{\link{results}} for more information about filtering using
#' Cook's distance, and the 'Dealing with outliers' section of the vignette.
#' Unlike the behavior of \code{\link{replaceOutliers}}, here original counts are
#' kept in the matrix returned by \code{\link{counts}}, original Cook's
#' distances are kept in \code{assays(dds)[["cooks"]]}, and the replacement
#' counts used for fitting are kept in \code{assays(dds)[["replaceCounts"]]}.
#'
#' Note that if a log2 fold change prior is used (betaPrior=TRUE)
#' then expanded model matrices will be used in fitting. These are
#' described in \code{\link{nbinomWaldTest}} and in the vignette. The
#' \code{contrast} argument of \code{\link{results}} should be used for
#' generating results tables.
#' 
#' @return a \code{\link{DESeqDataSet}} object with results stored as
#' metadata columns. These results should accessed by calling the \code{\link{results}}
#' function. By default this will return the log2 fold changes and p-values for the last
#' variable in the design formula.  See \code{\link{results}} for how to access results
#' for other variables.
#'
#' @param object a DESeqDataSet object, see the constructor functions
#' \code{\link{DESeqDataSet}},
#' \code{\link{DESeqDataSetFromMatrix}},
#' \code{\link{DESeqDataSetFromHTSeqCount}}.
#' @param test either "Wald" or "LRT", which will then use either 
#' Wald significance tests (defined by \code{\link{nbinomWaldTest}}),
#' or the likelihood ratio test on the difference in deviance between a
#' full and reduced model formula (defined by \code{\link{nbinomLRT}})
#' @param fitType either "parametric", "local", or "mean"
#' for the type of fitting of dispersions to the mean intensity.
#' See \code{\link{estimateDispersions}} for description.
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients 
#' See \code{\link{nbinomWaldTest}} for description of the calculation
#' of the beta prior. By default, the beta prior is used only for the
#' Wald test, but can also be specified for the likelihood ratio test.
#' @param full for \code{test="LRT"}, the full model formula,
#' which is restricted to the formula in \code{design(object)}.
#' alternatively, it can be a model matrix constructed by the user.
#' advanced use: specifying a model matrix for full and \code{test="Wald"}
#' is possible if \code{betaPrior=FALSE}
#' @param reduced for \code{test="LRT"}, a reduced formula to compare against,
#' i.e., the full formula with the term(s) of interest removed.
#' alternatively, it can be a model matrix constructed by the user
#' @param quiet whether to print messages at each step
#' @param minReplicatesForReplace the minimum number of replicates required
#' in order to use \code{\link{replaceOutliers}} on a
#' sample. If there are samples with so many replicates, the model will
#' be refit after these replacing outliers, flagged by Cook's distance.
#' Set to \code{Inf} in order to never replace outliers.
#' @param modelMatrixType either "standard" or "expanded", which describe
#' how the model matrix, X of the GLM formula is formed.
#' "standard" is as created by \code{model.matrix} using the
#' design formula. "expanded" includes an indicator variable for each
#' level of factors in addition to an intercept. for more information
#' see the Description of \code{\link{nbinomWaldTest}}.
#' betaPrior must be set to TRUE in order for expanded model matrices
#' to be fit.
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see next argument \code{BPPARAM}.
#' A note on running in parallel using \code{BiocParallel}: it may be
#' advantageous to remove large, unneeded objects from your current
#' R environment before calling \code{DESeq},
#' as it is possible that R's internal garbage collection
#' will copy these files while running on worker nodes.
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}.
#' If not specified, the parameters last registered with
#' \code{\link{register}} will be used.
#' 
#' @author Michael Love
#' 
#' @references
#'
#' Michael I Love, Wolfgang Huber, Simon Anders: Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 2014, 15:550. \url{http://dx.doi.org/10.1186/s13059-014-0550-8}
#' @import BiocGenerics BiocParallel S4Vectors IRanges GenomicRanges SummarizedExperiment Biobase Rcpp methods
#'
#' @importFrom locfit locfit
#' @importFrom genefilter rowVars filtered_p
#' @importFrom Hmisc wtd.quantile
#' 
#' @importFrom graphics axis hist plot points
#' @importFrom stats Gamma as.formula coefficients df dnbinom dnorm formula glm loess lowess model.matrix optim p.adjust pchisq pnorm prcomp predict pt qf qnorm rchisq relevel rnbinom rnorm runif splinefun terms terms.formula
#' @importFrom utils read.table
#' 
#' @useDynLib DESeq2
#'
#' @seealso \code{\link{nbinomWaldTest}}, \code{\link{nbinomLRT}}
#'
#' @examples
#'
#' # see vignette for suggestions on generating
#' # count tables from RNA-Seq data
#' cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
#' cond <- factor(rep(1:2, each=5))
#'
#' # object construction
#' dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
#'
#' # standard analysis
#' dds <- DESeq(dds)
#' res <- results(dds)
#'
#' # an alternate analysis: likelihood ratio test
#' ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
#' resLRT <- results(ddsLRT)
#'
#' @export
DESeq <- function(object, test=c("Wald","LRT"),
                  fitType=c("parametric","local","mean"), betaPrior,
                  full=design(object), reduced, quiet=FALSE,
                  minReplicatesForReplace=7, modelMatrixType,
                  parallel=FALSE, BPPARAM=bpparam()) {
  # check arguments
  stopifnot(is(object, "DESeqDataSet"))
  test <- match.arg(test, choices=c("Wald","LRT"))
  fitType <- match.arg(fitType, choices=c("parametric","local","mean"))
  stopifnot(is.logical(quiet))
  stopifnot(is.numeric(minReplicatesForReplace))
  stopifnot(is.logical(parallel))
  modelAsFormula <- !is.matrix(full)
  
  if (missing(betaPrior)) {
    betaPrior <- if (modelAsFormula) {
      termsOrder <- attr(terms.formula(design(object)),"order")
      interactionPresent <- any(termsOrder > 1)
      # use beta prior for Wald tests and when no interaction terms are included
      (test == "Wald") & !interactionPresent
    } else {
      FALSE
    }
  } else {
    stopifnot(is.logical(betaPrior))
  }
  # get rid of any NA in the mcols(mcols(object))
  object <- sanitizeRowRanges(object)
  
  if (test == "LRT") {
    if (missing(reduced)) {
      stop("likelihood ratio test requires a 'reduced' design, see ?DESeq")
    }
    if (!missing(modelMatrixType) && modelMatrixType == "expanded") {
      stop("test='LRT' only implemented for standard model matrices")
    }
    if (is.matrix(full) | is.matrix(reduced)) {
      if (!(is.matrix(full) & is.matrix(reduced))) {
        stop("if one of 'full' and 'reduced' is a matrix, the other must be also a matrix")
      }
    }
    if (modelAsFormula) {
      checkLRT(full, reduced)
    } else {
      checkFullRank(full)
      checkFullRank(reduced)
      if (ncol(full) <= ncol(reduced)) {
        stop("the number of columns of 'full' should be more than the number of columns of 'reduced'")
      }
    }
  }
  if (test == "Wald" & !missing(reduced)) {
    stop("'reduced' ignored when test='Wald'")
  }
  
  if (modelAsFormula) {

    # run some tests common to DESeq, nbinomWaldTest, nbinomLRT
    designAndArgChecker(object, betaPrior)
    
    if (full != design(object)) {
      stop("'full' specified as formula should equal design(object)")
    }
    modelMatrix <- NULL
  } else {
    if (betaPrior == TRUE) {
      stop("betaPrior=TRUE is not supported for user-provided model matrices")
    }
    checkFullRank(full)
    # this will be used for dispersion estimation and testing
    modelMatrix <- full
  }
 
  attr(object, "betaPrior") <- betaPrior
  stopifnot(length(parallel) == 1 & is.logical(parallel))
  
  if (!is.null(sizeFactors(object)) || !is.null(normalizationFactors(object))) {
    if (!quiet) {
      if (!is.null(normalizationFactors(object))) {
        message("using pre-existing normalization factors")
      } else {
        message("using pre-existing size factors")
      }
    }
  } else {
    if (!quiet) message("estimating size factors")
    object <- estimateSizeFactors(object)
  }
  
  if (!parallel) {
    if (!quiet) message("estimating dispersions")
    object <- estimateDispersions(object, fitType=fitType, quiet=quiet, modelMatrix=modelMatrix)
    if (!quiet) message("fitting model and testing")
    if (test == "Wald") {
      object <- nbinomWaldTest(object, betaPrior=betaPrior, quiet=quiet,
                               modelMatrix=modelMatrix,
                               modelMatrixType=modelMatrixType)
    } else if (test == "LRT") {
      object <- nbinomLRT(object, full=full, reduced=reduced,
                          betaPrior=betaPrior, quiet=quiet)
    }
  } else if (parallel) {
    object <- DESeqParallel(object, test=test, fitType=fitType,
                            betaPrior=betaPrior, full=full, reduced=reduced,
                            quiet=quiet, modelMatrix=modelMatrix,
                            modelMatrixType=modelMatrixType,
                            BPPARAM=BPPARAM)
  }

  # if there are sufficient replicates, then pass through to refitting function
  sufficientReps <- any(nOrMoreInCell(attr(object,"modelMatrix"),minReplicatesForReplace))
  if (sufficientReps) {
    object <- refitWithoutOutliers(object, test=test, betaPrior=betaPrior,
                                   full=full, reduced=reduced, quiet=quiet,
                                   minReplicatesForReplace=minReplicatesForReplace,
                                   modelMatrix=modelMatrix,
                                   modelMatrixType=modelMatrixType)
  }

  # stash the package version (again, also in construction)
  metadata(object)[["version"]] <- packageVersion("DESeq2")
  
  object
}

#' Make a simulated DESeqDataSet
#'
#' Constructs a simulated dataset of Negative Binomial data from
#' two conditions. By default, there are no fold changes between
#' the two conditions, but this can be adjusted with the \code{betaSD} argument.
#'
#' @param n number of rows
#' @param m number of columns
#' @param betaSD the standard deviation for non-intercept betas, i.e. beta ~ N(0,betaSD)
#' @param interceptMean the mean of the intercept betas (log2 scale)
#' @param interceptSD the standard deviation of the intercept betas (log2 scale)
#' @param dispMeanRel a function specifying the relationship of the dispersions on
#' \code{2^trueIntercept}
#' @param sizeFactors multiplicative factors for each sample
#'
#' @return a \code{\link{DESeqDataSet}} with true dispersion,
#' intercept and beta values in the metadata columns.  Note that the true
#' betas are provided on the log2 scale.
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' dds
#'
#' @export
makeExampleDESeqDataSet <- function(n=1000,m=12,betaSD=0,interceptMean=4,interceptSD=2,
                                    dispMeanRel=function(x) 4/x + .1,sizeFactors=rep(1,m)) {
  beta <- cbind(rnorm(n,interceptMean,interceptSD),rnorm(n,0,betaSD))
  dispersion <- dispMeanRel(2^(beta[,1]))
  colData <- DataFrame(condition=factor(rep(c("A","B"),times=c(ceiling(m/2),floor(m/2)))))
  x <- if (m > 1) {
    model.matrix(~ colData$condition)
  } else {
    cbind(rep(1,m),rep(0,m))
  }
  mu <- t(2^(x %*% t(beta)) * sizeFactors)
  countData <- matrix(rnbinom(m*n, mu=mu, size=1/dispersion), ncol=m)
  mode(countData) <- "integer"
  colnames(countData) <- paste("sample",1:m,sep="")
  rowRanges <- GRanges("1",IRanges(start=(1:n - 1) * 100 + 1,width=100))
  names(rowRanges) <- paste0("gene",1:n)

  # set environment to global environment,
  # to avoid the formula carrying with it all the objects
  # here including 'object' itself.
  design <- if (m > 1) {
    as.formula("~ condition", env=.GlobalEnv)
  } else {
    as.formula("~ 1", env=.GlobalEnv)
  }
  
  object <- DESeqDataSetFromMatrix(countData = countData,
                                   colData = colData,
                                   design = design,
                                   rowRanges = rowRanges)
  trueVals <- DataFrame(trueIntercept = beta[,1],
                        trueBeta = beta[,2],
                        trueDisp = dispersion)
  mcols(trueVals) <- DataFrame(type=rep("input",ncol(trueVals)),
                               description=c("simulated intercept values",
                                 "simulated beta values",
                                 "simulated dispersion values"))
  mcols(object) <- cbind(mcols(object),trueVals)
  return(object)
}


#' Low-level function to estimate size factors with robust regression.
#' 
#' Given a matrix or data frame of count data, this function estimates the size
#' factors as follows: Each column is divided by the geometric means of the
#' rows. The median (or, if requested, another location estimator) of these
#' ratios (skipping the genes with a geometric mean of zero) is used as the size
#' factor for this column. Typically, one will not call this function directly, but use
#' \code{\link{estimateSizeFactors}}.
#' 
#' @param counts a matrix or data frame of counts, i.e., non-negative integer
#' values
#' @param locfunc a function to compute a location for a sample. By default, the
#' median is used. However, especially for low counts, the
#' \code{\link[genefilter]{shorth}} function from genefilter may give better results.
#' @param geoMeans by default this is not provided, and the
#' geometric means of the counts are calculated within the function.
#' A vector of geometric means from another count matrix can be provided
#' for a "frozen" size factor calculation
#' @param controlGenes optional, numeric or logical index vector specifying those genes to
#' use for size factor estimation (e.g. housekeeping or spike-in genes)
#' @return a vector with the estimates size factors, one element per column
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactors}}
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet()
#' estimateSizeFactorsForMatrix(counts(dds))
#' geoMeans <- exp(rowMeans(log(counts(dds))))
#' estimateSizeFactorsForMatrix(counts(dds),geoMeans=geoMeans)
#' 
#' @export
estimateSizeFactorsForMatrix <- function( counts, locfunc = stats::median, geoMeans, controlGenes )
{
  if (missing(geoMeans)) {
    loggeomeans <- rowMeans(log(counts))
  } else {
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if ( !( is.numeric(controlGenes) | is.logical(controlGenes) ) ) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes,,drop=FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  sf
}

#' Low-level functions to fit dispersion estimates
#'
#' Normal users should instead use \code{\link{estimateDispersions}}.
#' These low-level functions are called by \code{\link{estimateDispersions}},
#' but are exported and documented for non-standard usage.
#' For instance, it is possible to replace fitted values with a custom fit and continue
#' with the maximum a posteriori dispersion estimation, as demonstrated in the
#' examples below.
#'
#' @param object a DESeqDataSet
#' @param fitType either "parametric", "local", or "mean"
#' for the type of fitting of dispersions to the mean intensity.
#' See \code{\link{estimateDispersions}} for description.
#' @param outlierSD the number of standard deviations of log
#' gene-wise estimates above the prior mean (fitted value),
#' above which dispersion estimates will be labelled
#' outliers. Outliers will keep their original value and
#' not be shrunk using the prior.
#' @param dispPriorVar the variance of the normal prior on the log dispersions.
#' If not supplied, this is calculated as the difference between
#' the mean squared residuals of gene-wise estimates to the
#' fitted dispersion and the expected sampling variance
#' of the log dispersion
#' @param minDisp small value for the minimum dispersion, to allow
#' for calculations in log scale, one order of magnitude above this value is used
#' as a test for inclusion in mean-dispersion fitting
#' @param kappa_0 control parameter used in setting the initial proposal
#' in backtracking search, higher kappa_0 results in larger steps
#' @param dispTol control parameter to test for convergence of log dispersion,
#' stop when increase in log posterior is less than dispTol
#' @param maxit control parameter: maximum number of iterations to allow for convergence
#' @param quiet whether to print messages at each step
#' @param modelMatrix for advanced use only,
#' a substitute model matrix for gene-wise and MAP dispersion estimation
#' @param niter number of times to iterate between estimation of means and
#' estimation of dispersion
#'
#' @return a DESeqDataSet with gene-wise, fitted, or final MAP
#' dispersion estimates in the metadata columns of the object.
#' 
#' \code{estimateDispersionsPriorVar} is called inside of \code{estimateDispersionsMAP}
#' and stores the dispersion prior variance as an attribute of
#' \code{dispersionFunction(dds)}, which can be manually provided to
#' \code{estimateDispersionsMAP} for parallel execution.
#'
#' @aliases estimateDispersionsGeneEst estimateDispersionsFit estimateDispersionsMAP estimateDispersionsPriorVar
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersionsGeneEst(dds)
#' dds <- estimateDispersionsFit(dds)
#' dds <- estimateDispersionsMAP(dds)
#' plotDispEsts(dds) 
#'
#' # after having run estimateDispersionsFit()
#' # the dispersion prior variance over all genes
#' # can be obtained like so:
#' 
#' dispPriorVar <- estimateDispersionsPriorVar(dds)
#' 
#' @seealso \code{\link{estimateDispersions}}
#'
#' @export
estimateDispersionsGeneEst <- function(object, minDisp=1e-8, kappa_0=1,
                                       dispTol=1e-6, maxit=100, quiet=FALSE,
                                       modelMatrix=NULL, niter=1) {
  if (!is.null(mcols(object)$dispGeneEst)) {
    if (!quiet) message("found already estimated gene-wise dispersions, removing these")
    removeCols <- c("dispGeneEst")
    mcols(object) <- mcols(object)[,!names(mcols(object)) %in% removeCols,drop=FALSE]
  }
  stopifnot(length(minDisp) == 1)
  stopifnot(length(kappa_0) == 1)
  stopifnot(length(dispTol) == 1)
  stopifnot(length(maxit) == 1)
  if (log(minDisp/10) <= -30) {
    stop("for computational stability, log(minDisp/10) should be above -30")
  }

  # in case the class of the mcols(mcols(object)) are not character
  object <- sanitizeRowRanges(object)

  if (is.null(modelMatrix)) {
    modelMatrix <- model.matrix(design(object), data=colData(object))
    checkFullRank(modelMatrix)
    if (nrow(modelMatrix) == ncol(modelMatrix)) {
      stop("the number of samples and the number of model coefficients are equal,
  i.e., there are no replicates to estimate the dispersion.
  use an alternate design formula")
    }
  } else {
    message("using supplied model matrix")
  }
  
  object <- getBaseMeansAndVariances(object)

  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  
  # this rough dispersion estimate (alpha_hat)
  # is for estimating mu
  # and for the initial starting point for line search
  # first check if model matrix is full rank
  fullRank <- qr(modelMatrix)$rank == ncol(modelMatrix)
  alpha_hat <- if (fullRank) {
    # if full rank use this estimator which compares normalized counts to mu
    roughDisp <- roughDispEstimate(y = counts(objectNZ,normalized=TRUE),
                                   x = modelMatrix)
    momentsDisp <- momentsDispEstimate(objectNZ)
    pmin(roughDisp, momentsDisp)
  } else {
    # if not full rank use method of moments across all samples
    momentsDispEstimate(objectNZ)
  }

  # bound the rough estimated alpha between minDisp and maxDisp for numeric stability
  maxDisp <- max(10, ncol(object))
  alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, alpha_hat), maxDisp)

  stopifnot(length(niter) == 1 & niter > 0)

  # below, iterate between mean and dispersion estimation (niter) times
  fitidx <- rep(TRUE,nrow(objectNZ))
  mu <- matrix(0, nrow=nrow(objectNZ), ncol=ncol(objectNZ))
  dispIter <- numeric(nrow(objectNZ))
  # bound the estimated count
  # this helps make the fitting more robust,
  # because 1/mu occurs in the weights for the NB GLM
  minmu <- 0.5
  for (iter in seq_len(niter)) {
    fit <- fitNbinomGLMs(objectNZ[fitidx,,drop=FALSE],
                         alpha_hat=alpha_hat[fitidx],
                         modelMatrix=modelMatrix)
    fitMu <- fit$mu
    fitMu[fitMu < minmu] <- minmu
    mu[fitidx,] <- fitMu
    # use of kappa_0 in backtracking search
    # initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
    # use log(minDisp/10) to stop if dispersions going to -infinity
    dispRes <- fitDispWrapper(ySEXP = counts(objectNZ)[fitidx,,drop=FALSE],
                              xSEXP = fit$modelMatrix,
                              mu_hatSEXP = fitMu,
                              log_alphaSEXP = log(alpha_hat)[fitidx],
                              log_alpha_prior_meanSEXP = log(alpha_hat)[fitidx],
                              log_alpha_prior_sigmasqSEXP = 1, min_log_alphaSEXP = log(minDisp/10),
                              kappa_0SEXP = kappa_0, tolSEXP = dispTol,
                              maxitSEXP = maxit, use_priorSEXP = FALSE)
    dispIter[fitidx] <- dispRes$iter
    alpha_hat_new[fitidx] <- pmin(exp(dispRes$log_alpha), maxDisp)
    # only rerun those rows which moved
    fitidx <- abs(log(alpha_hat_new) - log(alpha_hat)) > .05
    alpha_hat <- alpha_hat_new
    if (sum(fitidx) == 0) break
  }
 
  # dont accept moves if the log posterior did not
  # increase by more than one millionth,
  # and set the small estimates to the minimum dispersion
  dispGeneEst <- alpha_hat
  if (niter == 1) {
    noIncrease <- dispRes$last_lp < dispRes$initial_lp + abs(dispRes$initial_lp)/1e6
    dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
  }
  dispGeneEstConv <- dispIter < maxit
 
  # when lacking convergence from the C++ routine
  # we use an R function to estimate dispersions
  # by evaluating a grid of posterior evaluations
  refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp*10
  if (sum(refitDisp) > 0) {
    dispGrid <- fitDispGridWrapper(y=counts(objectNZ)[refitDisp,,drop=FALSE],x=modelMatrix,
                                   mu=mu[refitDisp,,drop=FALSE],
                                   logAlphaPriorMean=rep(0,sum(refitDisp)),
                                   logAlphaPriorSigmaSq=1,usePrior=FALSE)
    dispGeneEst[refitDisp] <- dispGrid
  }
  dispGeneEst <- pmin(pmax(dispGeneEst, minDisp), maxDisp)
  
  dispDataFrame <- buildDataFrameWithNARows(list(dispGeneEst=dispGeneEst),
                                            mcols(object)$allZero)
  mcols(dispDataFrame) <- DataFrame(type=rep("intermediate",ncol(dispDataFrame)),
                                    description=c("gene-wise estimates of dispersion"))
  mcols(object) <- cbind(mcols(object), dispDataFrame)
  assays(object)[["mu"]] <- buildMatrixWithNARows(mu, mcols(object)$allZero)
  
  return(object)
}

#' @rdname estimateDispersionsGeneEst
#' @export
estimateDispersionsFit <- function(object,fitType=c("parametric","local","mean"),
                                   minDisp=1e-8, quiet=FALSE) {

  if (is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  useForFit <- mcols(objectNZ)$dispGeneEst > 100*minDisp
  if (sum(useForFit) == 0) {
    stop("all gene-wise dispersion estimates are within 2 orders of magnitude
  from the minimum value, and so the standard curve fitting techniques will not work.
  One can instead use the gene-wise estimates as final estimates:
  dds <- estimateDispersionsGeneEst(dds)
  dispersions(dds) <- mcols(dds)$dispGeneEst
  ...then continue with testing using nbinomWaldTest or nbinomLRT")
  }
  
  fitType <- match.arg(fitType, choices=c("parametric","local","mean"))
  stopifnot(length(fitType)==1)
  stopifnot(length(minDisp)==1)
  if (fitType == "parametric") {
    trial <- try(dispFunction <- parametricDispersionFit(mcols(objectNZ)$baseMean[useForFit],
                                                         mcols(objectNZ)$dispGeneEst[useForFit]),
                 silent=TRUE)
    if (inherits(trial,"try-error")) {
      message("-- note: fitType='parametric', but the dispersion trend was not well captured by the
   function: y = a/x + b, and a local regression fit was automatically substituted.
   specify fitType='local' or 'mean' to avoid this message next time.")
      fitType <- "local"
    }
  }
  if (fitType == "local") {
    dispFunction <- localDispersionFit(means = mcols(objectNZ)$baseMean[useForFit],
                                       disps = mcols(objectNZ)$dispGeneEst[useForFit],
                                       minDisp = minDisp)
  }
  if (fitType == "mean") {
    useForMean <- mcols(objectNZ)$dispGeneEst > 10*minDisp
    meanDisp <- mean(mcols(objectNZ)$dispGeneEst[useForMean],na.rm=TRUE,trim=0.001)
    dispFunction <- function(means) meanDisp
    attr( dispFunction, "mean" ) <- meanDisp
  }
  if (!(fitType %in% c("parametric","local","mean"))) {
    stop("unknown fitType")
  }
 
  # store the dispersion function and attributes
  attr( dispFunction, "fitType" ) <- fitType
  if (quiet) {
    suppressMessages({ dispersionFunction(object) <- dispFunction })
  } else {
    dispersionFunction(object) <- dispFunction
  }
  
  return(object)
}

#' @rdname estimateDispersionsGeneEst
#' @export
estimateDispersionsMAP <- function(object, outlierSD=2, dispPriorVar,
                                   minDisp=1e-8, kappa_0=1, dispTol=1e-6,
                                   maxit=100, modelMatrix=NULL, quiet=FALSE) {
  stopifnot(length(outlierSD)==1)
  stopifnot(length(minDisp)==1)
  stopifnot(length(kappa_0)==1)
  stopifnot(length(dispTol)==1)
  stopifnot(length(maxit)==1)
  if (is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }
  if (!is.null(mcols(object)$dispersion)) {
    if (!quiet) message("found already estimated dispersions, removing these")
    removeCols <- c("dispersion","dispOutlier","dispMAP","dispIter","dispConv")
    mcols(object) <- mcols(object)[,!names(mcols(object)) %in% removeCols,drop=FALSE]
  }

  if (is.null(modelMatrix)) {
    modelMatrix <- model.matrix(design(object), data=colData(object))
  } else {
    message("using supplied model matrix")
  }
  
  # fill in the calculated dispersion prior variance
  if (missing(dispPriorVar)) {
    # if no gene-wise estimates above minimum
    if (sum(mcols(object)$dispGeneEst >= minDisp*100,na.rm=TRUE) == 0) {
      warning(paste0("all genes have dispersion estimates < ",minDisp*10,
                     ", returning disp = ",minDisp*10))
      resultsList <- list(dispersion = rep(minDisp*10, sum(!mcols(object)$allZero)))
      dispDataFrame <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
      mcols(dispDataFrame) <- DataFrame(type="intermediate",
                                        description="final estimates of dispersion")
      mcols(object) <- cbind(mcols(object), dispDataFrame)
      dispFn <- dispersionFunction(object)
      attr( dispFn, "dispPriorVar" ) <- 0.25
      dispersionFunction(object, estimateVar=FALSE) <- dispFn
      return(object)
    }
    dispPriorVar <- estimateDispersionsPriorVar(object, modelMatrix=modelMatrix)
    dispFn <- dispersionFunction(object)
    attr( dispFn, "dispPriorVar" ) <- dispPriorVar
    dispersionFunction(object, estimateVar=FALSE) <- dispFn
  } else {
    dispFn <- dispersionFunction(object)
    attr( dispFn, "dispPriorVar" ) <- dispPriorVar
    dispersionFunction(object, estimateVar=FALSE) <- dispFn
  }

  stopifnot(length(dispPriorVar)==1)

  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  varLogDispEsts <- attr( dispersionFunction(object), "varLogDispEsts" )
  
  # set prior variance for fitting dispersion
  log_alpha_prior_sigmasq <- dispPriorVar

  # get previously calculated mu
  mu <- assays(objectNZ)[["mu"]]

  # start fitting at gene estimate unless the points are one order of magnitude
  # below the fitted line, then start at fitted line
  dispInit <- ifelse(mcols(objectNZ)$dispGeneEst >  0.1 * mcols(objectNZ)$dispFit,
                     mcols(objectNZ)$dispGeneEst,
                     mcols(objectNZ)$dispFit)

  # if any missing values, fill in the fitted value to initialize
  dispInit[is.na(dispInit)] <- mcols(objectNZ)$dispFit[is.na(dispInit)]
  
  # run with prior
  dispResMAP <- fitDispWrapper(ySEXP = counts(objectNZ), xSEXP = modelMatrix, mu_hatSEXP = mu,
                               log_alphaSEXP = log(dispInit),
                               log_alpha_prior_meanSEXP = log(mcols(objectNZ)$dispFit),
                               log_alpha_prior_sigmasqSEXP = log_alpha_prior_sigmasq,
                               min_log_alphaSEXP = log(minDisp/10),
                               kappa_0SEXP = kappa_0, tolSEXP = dispTol,
                               maxitSEXP = maxit, use_priorSEXP = TRUE)

  # prepare dispersions for storage in mcols(object)
  dispMAP <- exp(dispResMAP$log_alpha) 

  # when lacking convergence from the C++ routine
  # we use an R function to estimate dispersions.
  # This finds the maximum of a smooth curve along a
  # grid of posterior evaluations
  dispConv <- dispResMAP$iter < maxit
  refitDisp <- !dispConv
  if (sum(refitDisp) > 0) {
    dispInR <- fitDispGridWrapper(y = counts(objectNZ)[refitDisp,,drop=FALSE], x = modelMatrix,
                                  mu = mu[refitDisp,,drop=FALSE],
                                  logAlphaPriorMean = log(mcols(objectNZ)$dispFit)[refitDisp],
                                  logAlphaPriorSigmaSq = log_alpha_prior_sigmasq,
                                  usePrior=TRUE)
    dispMAP[refitDisp] <- dispInR
  }

  # bound the dispersion estimate between minDisp and maxDisp for numeric stability
  maxDisp <- max(10, ncol(object))
  dispMAP <- pmin(pmax(dispMAP, minDisp), maxDisp)
  
  dispersionFinal <- dispMAP
  
  # detect outliers which have gene-wise estimates
  # outlierSD * standard deviation of log gene-wise estimates
  # above the fitted mean (prior mean)
  # and keep the original gene-est value for these.
  # Note: we use the variance of log dispersions estimates
  # from all the genes, not only those from below
  dispOutlier <- log(mcols(objectNZ)$dispGeneEst) >
                 log(mcols(objectNZ)$dispFit) +
                 outlierSD * sqrt(varLogDispEsts)
  dispOutlier[is.na(dispOutlier)] <- FALSE
  dispersionFinal[dispOutlier] <- mcols(objectNZ)$dispGeneEst[dispOutlier]
 
  resultsList <- list(dispersion = dispersionFinal,
                      dispIter = dispResMAP$iter,
                      dispOutlier = dispOutlier,
                      dispMAP = dispMAP)

  dispDataFrame <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
  mcols(dispDataFrame) <- DataFrame(type=rep("intermediate",ncol(dispDataFrame)),
                                    description=c("final estimate of dispersion",
                                      "number of iterations",
                                      "dispersion flagged as outlier",
                                      "maximum a posteriori estimate"))

  mcols(object) <- cbind(mcols(object), dispDataFrame)
  return(object)
}

#' @rdname estimateDispersionsGeneEst
#' @export
estimateDispersionsPriorVar <- function(object, minDisp=1e-8, modelMatrix=NULL) {
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  aboveMinDisp <- mcols(objectNZ)$dispGeneEst >= minDisp*100
  if (is.null(modelMatrix)) {
    modelMatrix <- model.matrix(design(object), data=colData(object))
  }
  # estimate the variance of the distribution of the
  # log dispersion estimates around the fitted value
  dispResiduals <- log(mcols(objectNZ)$dispGeneEst) - log(mcols(objectNZ)$dispFit)
  if (sum(aboveMinDisp,na.rm=TRUE) == 0) {
    stop("no data found which is greater than minDisp")
  }
  
  varLogDispEsts <- attr(dispersionFunction(object), "varLogDispEsts")
  
  m <- nrow(modelMatrix)
  p <- ncol(modelMatrix)

  # if the residual degrees of freedom is between 1 and 3, the distribution
  # of log dispersions is especially asymmetric and poorly estimated
  # by the MAD. we then use an alternate estimator, a monte carlo
  # approach to match the distribution
  if (((m - p) <= 3) & (m > p)) {
    # in order to produce identical results we set the seed, 
    # and so we need to save and restore the .Random.seed value first
    if (exists(".Random.seed")) {
      oldRandomSeed <- .Random.seed
    }
    set.seed(2)
    # The residuals are the observed distribution we try to match
    obsDist <- dispResiduals[aboveMinDisp]
    brks <- -20:20/2
    obsDist <- obsDist[obsDist > min(brks) & obsDist < max(brks)]
    obsVarGrid <- seq(from=0,to=8,length=200)
    obsDistHist <- hist(obsDist,breaks=brks,plot=FALSE)
    klDivs <- sapply(obsVarGrid, function(x) {
      randDist <- log(rchisq(1e4,df=(m-p))) + rnorm(1e4,0,sqrt(x)) - log(m - p)
      randDist <- randDist[randDist > min(brks) & randDist < max(brks)]
      randDistHist <- hist(randDist,breaks=brks,plot=FALSE)
      z <- c(obsDistHist$density,randDistHist$density)
      small <- min(z[z > 0])
      kl <- sum(obsDistHist$density * (log(obsDistHist$density + small) - log(randDistHist$density + small)))
      kl
    })
    lofit <- loess(klDivs ~ obsVarGrid, span=.2)
    obsVarFineGrid <- seq(from=0,to=8,length=1000)
    lofitFitted <- predict(lofit,obsVarFineGrid)
    argminKL <- obsVarFineGrid[which.min(lofitFitted)]
    expVarLogDisp <- trigamma((m - p)/2)
    dispPriorVar <- pmax(argminKL, 0.25)
    # finally, restore the .Random.seed if it existed beforehand
    if (exists("oldRandomSeed")) {
      .Random.seed <<- oldRandomSeed
    }

    return(dispPriorVar)
  }

  # estimate the expected sampling variance of the log estimates
  # Var(log(cX)) = Var(log(X))
  # X ~ chi-squared with m - p degrees of freedom
  if (m > p) {
    expVarLogDisp <- trigamma((m - p)/2)
    # set the variance of the prior using these two estimates
    # with a minimum of .25
    dispPriorVar <- pmax((varLogDispEsts - expVarLogDisp), 0.25)
  } else {
    # we have m = p, so do not try to subtract sampling variance
    dispPriorVar <- varLogDispEsts
    expVarLogDisp <- 0
  }

  dispPriorVar
}



#' Wald test for the GLM coefficients
#' 
#' This function tests for significance of coefficients in a Negative
#' Binomial GLM, using previously calculated \code{\link{sizeFactors}}
#' (or \code{\link{normalizationFactors}})
#' and dispersion estimates.  See \code{\link{DESeq}} for the GLM formula.
#' 
#' The fitting proceeds as follows: standard maximum likelihood estimates
#' for GLM coefficients (synonymous with "beta", "log2 fold change", "effect size")
#' are calculated. A zero-centered Normal prior distribution 
#' is assumed for the coefficients other than the intercept.
#' The variance of the prior distribution for each
#' non-intercept coefficient is calculated using the observed
#' distribution of the maximum likelihood coefficients.  
#' The final coefficients are then maximum a posteriori estimates
#' using this prior (Tikhonov/ridge regularization). See below for details on the
#' prior variance and the Methods section of the DESeq2 manuscript for more detail.
#' The use of a prior has little effect on genes with high counts and helps to
#' moderate the large spread in coefficients for genes with low counts.
#' For calculating Wald test p-values, the coefficients are scaled by their
#' standard errors and then compared to a standard Normal distribution. 
#'
#' The prior variance is calculated by matching the 0.05 upper quantile
#' of the observed MLE coefficients to a zero-centered Normal distribution.
#' In a change of methods since the 2014 paper,
#' the weighted upper quantile is calculated using the
#' \code{wtd.quantile} function from the Hmisc package. The weights are
#' the inverse of the expected variance of log counts, so the inverse of
#' \eqn{1/\bar{\mu} + \alpha_{tr}}{1/mu-bar + alpha_tr} using the mean of
#' normalized counts and the trended dispersion fit. The weighting ensures
#' that noisy estimates of log fold changes from small count genes do not
#' overly influence the calculation of the prior variance.
#' See \code{\link{estimateBetaPriorVar}}.
#' The final prior variance for a factor level is the average of the
#' estimated prior variance over all contrasts of all levels of the factor. 
#' Another change since the 2014 paper: when interaction terms are present
#' in the design, the prior on log fold changes is turned off
#' (for more details, see the vignette section, "Methods changes since
#' the 2014 DESeq2 paper").
#' 
#' When a log2 fold change prior is used (betaPrior=TRUE),
#' then \code{nbinomWaldTest} will by default use expanded model matrices,
#' as described in the \code{modelMatrixType} argument, unless this argument
#' is used to override the default behavior.
#' This ensures that log2 fold changes will be independent of the choice
#' of reference level. In this case, the beta prior variance for each factor
#' is calculated as the average of the mean squared maximum likelihood
#' estimates for each level and every possible contrast. The \code{\link{results}}
#' function without any arguments will automatically perform a contrast of the
#' last level of the last variable in the design formula over the first level.
#' The \code{contrast} argument of the \code{\link{results}} function can be used
#' to generate other comparisons.
#'  
#' The Wald test can be replaced with the \code{\link{nbinomLRT}}
#' for an alternative test of significance.
#'
#' @param object a DESeqDataSet
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients
#' @param betaPriorVar a vector with length equal to the number of
#' model terms including the intercept.
#' betaPriorVar gives the variance of the prior on the sample betas
#' on the log2 scale. if missing (default) this is estimated from the data
#' @param modelMatrix an optional matrix, typically this is set to NULL
#' and created within the function. only can be supplied if betaPrior=FALSE
#' @param modelMatrixType either "standard" or "expanded", which describe
#' how the model matrix, X of the formula in \code{\link{DESeq}}, is
#' formed. "standard" is as created by \code{model.matrix} using the
#' design formula. "expanded" includes an indicator variable for each
#' level of factors in addition to an intercept.
#' betaPrior must be set to TRUE in order for expanded model matrices
#' to be fit.
#' @param maxit the maximum number of iterations to allow for convergence of the
#' coefficient vector
#' @param useOptim whether to use the native optim function on rows which do not
#' converge within maxit
#' @param quiet whether to print messages at each step
#' @param useT whether to use a t-distribution as a null distribution,
#' for significance testing of the Wald statistics.
#' If FALSE, a standard normal null distribution is used.
#' @param df the degrees of freedom for the t-distribution
#' @param useQR whether to use the QR decomposition on the design
#' matrix X while fitting the GLM
#'
#' @return a DESeqDataSet with results columns accessible
#' with the \code{\link{results}} function.  The coefficients and standard errors are
#' reported on a log2 scale.
#'
#' @seealso \code{\link{DESeq}}, \code{\link{nbinomLRT}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' dds <- nbinomWaldTest(dds)
#' res <- results(dds)
#'
#' @export
nbinomWaldTest <- function(object, betaPrior, betaPriorVar,
                           modelMatrix=NULL, modelMatrixType,
                           maxit=100, useOptim=TRUE, quiet=FALSE,
                           useT=FALSE, df, useQR=TRUE) {
  if (is.null(dispersions(object))) {
    stop("testing requires dispersion estimates, first call estimateDispersions()")
  }
  stopifnot(length(maxit)==1)
  # in case the class of the mcols(mcols(object)) are not character
  object <- sanitizeRowRanges(object)
  
  if ("results" %in% mcols(mcols(object))$type) {
    if (!quiet) message("found results columns, replacing these")
    object <- removeResults(object)
  }
  if (is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }
  
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]

  if (is.null(modelMatrix)) {
    modelAsFormula <- TRUE
    termsOrder <- attr(terms.formula(design(object)),"order")
    interactionPresent <- any(termsOrder > 1)
    if (missing(betaPrior)) {
      betaPrior <- !interactionPresent
    }

    # run some tests common to DESeq, nbinomWaldTest, nbinomLRT
    designAndArgChecker(object, betaPrior)

    # what kind of model matrix to use
    stopifnot(is.logical(betaPrior))
    blindDesign <- design(object) == formula(~ 1)
    if (blindDesign) {
      betaPrior <- FALSE
    }
    if (missing(modelMatrixType) || is.null(modelMatrixType)) {
      modelMatrixType <- if (betaPrior) {
        "expanded"
      } else {
        "standard"
      }
    }
    if (modelMatrixType == "expanded" & !betaPrior) {
      stop("expanded model matrices require a beta prior")
    }
    # store modelMatrixType so it can be accessed by estimateBetaPriorVar
    attr(object, "modelMatrixType") <- modelMatrixType
    hasIntercept <- attr(terms(design(object)),"intercept") == 1
    renameCols <- hasIntercept
  } else {
    if (missing(betaPrior)) {
      betaPrior <- FALSE
    } else {
      if (betaPrior) stop("the model matrix can only be user-supplied if betaPrior=FALSE")
    } 
    message("using supplied model matrix")
    modelAsFormula <- FALSE
    attr(object, "modelMatrixType") <- "user-supplied"
    renameCols <- FALSE
  }

  if (!betaPrior) {
    # fit the negative binomial GLM without a prior
    # (in actuality a very wide prior with standard deviation 1e3 on log2 fold changes)
    fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim, useQR=useQR,
                         renameCols=renameCols, modelMatrix=modelMatrix)
    H <- fit$hat_diagonals
    modelMatrix <- fit$modelMatrix
    modelMatrixNames <- fit$modelMatrixNames
    # record the wide prior variance which was used in fitting
    betaPriorVar <- rep(1e6, ncol(fit$modelMatrix))
  } else {
    priorFitList <- fitGLMsWithPrior(object=object,
                                     maxit=maxit, useOptim=useOptim, useQR=useQR,
                                     betaPriorVar=betaPriorVar)
    fit <- priorFitList$fit
    H <- priorFitList$H
    betaPriorVar <- priorFitList$betaPriorVar
    modelMatrix <- priorFitList$modelMatrix
    mleBetaMatrix <- priorFitList$mleBetaMatrix

    # will add the MLE betas, so remove any which exist already
    # (possibly coming from estimateMLEForBetaPriorVar)
    mcols(object) <- mcols(object)[,grep("MLE_",names(mcols(object)),invert=TRUE)]
  }

  # store mu in case the user did not call estimateDispersionsGeneEst
  dimnames(fit$mu) <- NULL
  assays(objectNZ)[["mu"]] <- fit$mu
  assays(object)[["mu"]] <- buildMatrixWithNARows(fit$mu, mcols(object)$allZero)

  # store the prior variance directly as an attribute
  # of the DESeqDataSet object, so it can be pulled later by
  # the results function (necessary for setting max Cook's distance)
  attr(object,"betaPrior") <- betaPrior
  attr(object,"betaPriorVar") <- betaPriorVar
  attr(object,"modelMatrix") <- modelMatrix
  attr(object,"test") <- "Wald"

  # calculate Cook's distance
  dispModelMatrix <- if (modelAsFormula) {
    model.matrix(design(object), data=colData(object))
  } else {
    modelMatrix
  }
  attr(object,"dispModelMatrix") <- dispModelMatrix
  cooks <- calculateCooksDistance(objectNZ, H, dispModelMatrix)

  # record maximum Cook's
  maxCooks <- recordMaxCooks(design(object), colData(object), dispModelMatrix, cooks, nrow(objectNZ))

  # store Cook's distance for each sample
  assays(object)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)
  
  # add betas, standard errors and Wald p-values to the object
  modelMatrixNames <- colnames(modelMatrix)
  betaMatrix <- fit$betaMatrix
  colnames(betaMatrix) <- modelMatrixNames
  betaSE <- fit$betaSE
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)
  WaldStatistic <- betaMatrix/betaSE
  colnames(WaldStatistic) <- paste0("WaldStatistic_",modelMatrixNames)
  
  # if useT is set to TRUE, use a t-distribution
  if (useT) {
    dispPriorVar <- attr( dispersionFunction(object), "dispPriorVar" )
    stopifnot(length(df)==1)
    WaldPvalue <- 2*pt(abs(WaldStatistic),df=df,lower.tail=FALSE)
  } else {
    WaldPvalue <- 2*pnorm(abs(WaldStatistic),lower.tail=FALSE)
  }
  colnames(WaldPvalue) <- paste0("WaldPvalue_",modelMatrixNames)
  
  betaConv <- fit$betaConv

  if (any(!betaConv)) {
    if (!quiet) message(paste(sum(!betaConv),"rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"))
  }

  mleBetas <- if (betaPrior) {
    matrixToList(mleBetaMatrix)
  } else {
    NULL
  }
  
  resultsList <- c(matrixToList(betaMatrix),
                   matrixToList(betaSE),
                   mleBetas,
                   matrixToList(WaldStatistic),
                   matrixToList(WaldPvalue),
                   list(betaConv = betaConv,
                        betaIter = fit$betaIter,
                        deviance = -2 * fit$logLike,
                        maxCooks = maxCooks))
  
  WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
  
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)

  lfcType <- if (attr(object,"betaPrior")) "MAP" else "MLE"
  coefInfo <- paste(paste0("log2 fold change (",lfcType,"):"),modelMatrixNamesSpaces)
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)
  mleInfo <- if (betaPrior) {
    gsub("_"," ",colnames(mleBetaMatrix))
  } else {
    NULL
  }
  statInfo <- paste("Wald statistic:",modelMatrixNamesSpaces)
  pvalInfo <- paste("Wald test p-value:",modelMatrixNamesSpaces)

  mcols(WaldResults) <- DataFrame(type = rep("results",ncol(WaldResults)),
                                  description = c(coefInfo, seInfo, mleInfo, statInfo, pvalInfo,
                                    "convergence of betas",
                                    "iterations for betas",
                                    "deviance for the fitted model",
                                    "maximum Cook's distance for row"))
 
  mcols(object) <- cbind(mcols(object),WaldResults)
  return(object)
}



#' Steps for estimating the beta prior variance
#'
#' These lower-level functions are called within \code{\link{DESeq}} or \code{\link{nbinomWaldTest}}.
#' End users should use those higher-level function instead.
#' NOTE: \code{estimateBetaPriorVar} returns a numeric vector, not a DESEqDataSet!
#' For advanced users: to use these functions, first run \code{estimateMLEForBetaPriorVar}
#' and then run \code{estimateBetaPriorVar}.
#'
#' @param object a DESeqDataSet
#'
#' @param maxit as defined in \code{link{nbinomWaldTest}}
#' @param useOptim as defined in \code{link{nbinomWaldTest}}
#' @param useQR as defined in \code{link{nbinomWaldTest}}
#' 
#' @param betaPriorMethod the method for calculating the beta prior variance,
#' either "quanitle" or "weighted":
#' "quantile" matches a normal distribution using the upper quantile of the finite MLE betas.
#' "weighted" matches a normal distribution using the upper quantile, but weighting by the variance of the MLE betas.
#' @param upperQuantile the upper quantile to be used for the
#' "quantile" or "weighted" method of beta prior variance estimation
#'
#' 
#' @return for \code{estimateMLEForBetaPriorVar}, a DESeqDataSet, with the
#' necessary information stored in order to calculate the prior variance.
#' for \code{estimateBetaPriorVar}, the vector of variances for the prior
#' on the betas in the \code{\link{DESeq}} GLM
#'
#' @aliases estimateBetaPriorVar estimateMLEForBetaPriorVar
#' 
#' @export
estimateBetaPriorVar <- function(object, 
                                 betaPriorMethod=c("weighted","quantile"),
                                 upperQuantile=0.05) {
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]

  betaMatrix <- as.matrix(mcols(objectNZ)[,grep("MLE_", names(mcols(object))),drop=FALSE])
  colnamesBM <- colnames(betaMatrix)
  colnamesBM <- gsub("MLE_(.*)","\\1",colnamesBM)

  # renaming in reverse:
  # make these standard colnames as from model.matrix()
  convertNames <- renameModelMatrixColumns(colData(object),
                                           design(object))
  colnamesBM <- sapply(colnamesBM, function(x) {
    if (x %in% convertNames$to) {
      convertNames$from[convertNames$to == x]
    } else {
      x
    }
  })
  colnames(betaMatrix) <- colnamesBM
  
  # this is the model matrix from an MLE run
  modelMatrix <- model.matrix(design(objectNZ), colData(object))

  modelMatrixType <- attr(object, "modelMatrixType")
  
  betaPriorMethod <- match.arg(betaPriorMethod, choices=c("weighted","quantile"))

  # estimate the variance of the prior on betas
  # if expanded, first calculate LFC for all possible contrasts
  if (modelMatrixType == "expanded") {
    betaMatrix <- addAllContrasts(objectNZ, betaMatrix)
  }

  # weighting by 1/Var(log(K))
  # Var(log(K)) ~ Var(K)/mu^2 = 1/mu + alpha
  # and using the fitted alpha
  dispFit <- mcols(objectNZ)$dispFit
  if (is.null(dispFit)) {
    # betaPrior routine could have been called w/o the dispersion fitted trend
    dispFit <- mean(dispersions(objectNZ))
  }
  varlogk <- 1/mcols(objectNZ)$baseMean + dispFit
  weights <- 1/varlogk
  
  betaPriorVar <- if (nrow(betaMatrix) > 1) {
    apply(betaMatrix, 2, function(x) {
      # this test removes genes which have betas
      # tending to +/- infinity
      useFinite <- abs(x) < 10
      # if no more betas pass test, return wide prior
      if (sum(useFinite) == 0 ) {
        return(1e6)
      } else {
        if (betaPriorMethod=="quantile") {
          return(matchUpperQuantileForVariance(x[useFinite],upperQuantile))
        } else if (betaPriorMethod=="weighted") {
          return(matchWeightedUpperQuantileForVariance(x[useFinite],weights[useFinite],upperQuantile))
        }
      }
    })
  } else {
    (betaMatrix)^2
  }
  names(betaPriorVar) <- colnames(betaMatrix)

  # pre-v1.10 code for interactions and beta prior:
  # ------------------------------------------------------
  # find the names of betaPriorVar which correspond
  # to non-interaction terms and set these to a wide prior
  ## termsOrder <- attr(terms.formula(design(object)),"order")
  ## interactionPresent <- any(termsOrder > 1)  
  ## if (interactionPresent) {
  ##   nonInteractionCols <- getNonInteractionColumnIndices(objectNZ, modelMatrix)
  ##   if (modelMatrixType == "standard") widePrior <- 1e6 else widePrior <- 1e3
  ##   betaPriorVar[nonInteractionCols] <- widePrior
  ##   if (modelMatrixType == "expanded") {
  ##     # also set a wide prior for additional contrasts which were added
  ##     # for calculation of the prior variance in the case of
  ##     # expanded model matrices
  ##     designFactors <- getDesignFactors(objectNZ)
  ##     betaPriorVar[which(names(betaPriorVar) %in% paste0(designFactors,"Cntrst"))] <- widePrior
  ##   }
  ## }
  
  # intercept set to wide prior
  if ("Intercept" %in% names(betaPriorVar)) {
    betaPriorVar[which(names(betaPriorVar) == "Intercept")] <- 1e6
  }
  
  if (modelMatrixType == "expanded") {
    # bring over beta priors from the GLM fit without prior.
    # for factors: prior variance of each level are the average of the
    # prior variances for the levels present in the previous GLM fit
    betaPriorExpanded <- averagePriorsOverLevels(objectNZ, betaPriorVar)
    betaPriorVar <- betaPriorExpanded
  }
  
  betaPriorVar
}

#' @rdname estimateBetaPriorVar
#' @export
estimateMLEForBetaPriorVar <- function(object, maxit=100, useOptim=TRUE, useQR=TRUE) {
  # this function copies code from other functions,
  # in order to allow parallelization  
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]

  # this code copied from nbinomWaldTest()
  blindDesign <- design(object) == formula(~ 1)
  mmTypeTest <- !blindDesign
  modelMatrixType <- if (mmTypeTest) {
    "expanded"
  } else {
    "standard"
  }
  attr(object, "modelMatrixType") <- modelMatrixType

  # this code copied from fitGLMsWithPrior()
  fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim, useQR=useQR,
                       renameCols = (modelMatrixType == "standard"))
  modelMatrix <- fit$modelMatrix
  modelMatrixNames <- colnames(modelMatrix)
  H <- fit$hat_diagonal
  betaMatrix <- fit$betaMatrix
 
  modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  modelMatrixNames <- make.names(modelMatrixNames)
  colnames(betaMatrix) <- modelMatrixNames
  
  convertNames <- renameModelMatrixColumns(colData(object),
                                           design(objectNZ))
  convertNames <- convertNames[convertNames$from %in% modelMatrixNames,,drop=FALSE]
  modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
  mleBetaMatrix <- fit$betaMatrix
  colnames(mleBetaMatrix) <- paste0("MLE_",modelMatrixNames)
  # remove any MLE columns if they exist
  mcols(object) <- mcols(object)[,grep("MLE_",names(mcols(object)),invert=TRUE)]
  mcols(object) <- cbind(mcols(object), buildDataFrameWithNARows(DataFrame(mleBetaMatrix), mcols(object)$allZero))
  assays(object)[["H"]] <- buildMatrixWithNARows(H, mcols(object)$allZero)
  object
}




#' Likelihood ratio test (chi-squared test) for GLMs
#'
#' This function tests for significance of change in deviance between a
#' full and reduced model which are provided as \code{formula}.
#' Fitting uses previously calculated \code{\link{sizeFactors}} (or \code{\link{normalizationFactors}})
#' and dispersion estimates.
#' 
#' The difference in deviance is compared to a chi-squared distribution
#' with df = (reduced residual degrees of freedom - full residual degrees of freedom).
#' This function is comparable to the \code{nbinomGLMTest} of the previous version of DESeq
#' and an alternative to the default \code{\link{nbinomWaldTest}}.
#' 
#' @param object a DESeqDataSet
#' @param full the full model formula, this should be the formula in
#' \code{design(object)}.
#' alternatively, can be a matrix
#' @param reduced a reduced formula to compare against, e.g.
#' the full model with a term or terms of interest removed.
#' alternatively, can be a matrix
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients 
#' While the beta prior is used typically, for the Wald test, it can
#' also be specified for the likelihood ratio test. For more details
#' on the calculation, see \code{\link{nbinomWaldTest}}.
#' @param betaPriorVar a vector with length equal to the number of
#' model terms including the intercept.
#  betaPriorVar gives the variance of the prior on the sample betas,
#' which if missing is estimated from the rows which do not have any
#' zeros
#' @param maxit the maximum number of iterations to allow for convergence of the
#' coefficient vector
#' @param useOptim whether to use the native optim function on rows which do not
#' converge within maxit
#' @param quiet whether to print messages at each step
#' @param useQR whether to use the QR decomposition on the design
#' matrix X while fitting the GLM
#'
#' @return a DESeqDataSet with new results columns accessible
#' with the \code{\link{results}} function.  The coefficients and standard errors are
#' reported on a log2 scale.
#' 
#' @seealso \code{\link{DESeq}}, \code{\link{nbinomWaldTest}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet()
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' dds <- nbinomLRT(dds, reduced = ~ 1)
#' res <- results(dds)
#'
#' @export
nbinomLRT <- function(object, full=design(object), reduced,
                      betaPrior=FALSE, betaPriorVar,
                      maxit=100, useOptim=TRUE, quiet=FALSE,
                      useQR=TRUE) {

  if (is.null(dispersions(object))) {
    stop("testing requires dispersion estimates, first call estimateDispersions()")
  }
  if (missing(reduced)) {
    stop("provide a reduced formula for the LRT, e.g. nbinomLRT(object, reduced= ~1)")
  }

  # in case the class of the mcols(mcols(object)) are not character
  object <- sanitizeRowRanges(object)
  
  # run check on the formula
  modelAsFormula <- !(is.matrix(full) & is.matrix(reduced))
  if (modelAsFormula) {
    checkLRT(full, reduced)

    # run some tests common to DESeq, nbinomWaldTest, nbinomLRT
    designAndArgChecker(object, betaPrior)
    
    # try to form model matrices, test for difference
    # in residual degrees of freedom
    fullModelMatrix <- model.matrix(full,
                                    data=colData(object))
    reducedModelMatrix <- model.matrix(reduced,
                                       data=colData(object))
    df <- ncol(fullModelMatrix) - ncol(reducedModelMatrix)
  } else {
    if (betaPrior) {
      stop("user-supplied model matrices require betaPrior=FALSE")
    }
    message("using supplied model matrix")
    df <- ncol(full) - ncol(reduced)
  }
  
  if (df < 1) stop("less than one degree of freedom, perhaps full and reduced models are not in the correct order")
  
  if (any(mcols(mcols(object))$type == "results")) {
    if (!quiet) message("found results columns, replacing these")
    object <- removeResults(object)
  } 

  if (is.null(mcols(object)$allZero)) {
    object <- getBaseMeansAndVariances(object)
  }
  
  stopifnot(is.logical(betaPrior))
  if (modelAsFormula) {
    modelMatrixType <- "standard"
    # check for intercept
    hasIntercept <- attr(terms(design(object)),"intercept") == 1
    renameCols <- hasIntercept
  } else {
    modelMatrixType <- "user-supplied"
    renameCols <- FALSE
  }

  # store modelMatrixType so it can be accessed by estimateBetaPriorVar
  attr(object,"modelMatrixType") <- modelMatrixType

  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]

  if (!betaPrior) {
    if (modelAsFormula) {
      fullModel <- fitNbinomGLMs(objectNZ, modelFormula=full,
                                 renameCols=renameCols, maxit=maxit,
                                 useOptim=useOptim, useQR=useQR, warnNonposVar=FALSE)
      modelMatrix <- fullModel$modelMatrix
      reducedModel <- fitNbinomGLMs(objectNZ, modelFormula=reduced, maxit=maxit,
                                    useOptim=useOptim, useQR=useQR, warnNonposVar=FALSE)
    } else {
      fullModel <- fitNbinomGLMs(objectNZ, modelMatrix=full,
                                 renameCols=FALSE, maxit=maxit,
                                 useOptim=useOptim, useQR=useQR, warnNonposVar=FALSE)
      modelMatrix <- full
      reducedModel <- fitNbinomGLMs(objectNZ, modelMatrix=reduced,
                                    renameCols=FALSE, maxit=maxit,
                                    useOptim=useOptim, useQR=useQR, warnNonposVar=FALSE)
    }
    betaPriorVar <- rep(1e6, ncol(modelMatrix))
  } else {
    priorFull <- fitGLMsWithPrior(object=object,
                                  maxit=maxit, useOptim=useOptim, useQR=useQR,
                                  betaPriorVar=betaPriorVar)
    fullModel <- priorFull$fit
    modelMatrix <- fullModel$modelMatrix
    betaPriorVar <- priorFull$betaPriorVar
    mleBetaMatrix <- priorFull$mleBetaMatrix
    # form a reduced model matrix:
    # first find the dropped terms
    # then remove columns from the full model matrix which are
    # assigned to these terms
    fullModelTerms <- attr(terms(full),"term.labels")
    reducedModelTerms <- attr(terms(reduced),"term.labels")
    droppedTerms <- which(!fullModelTerms %in% reducedModelTerms)
    fullAssign <- attr(modelMatrix,"assign")
    idx <- !fullAssign %in% droppedTerms
    # now subsetting the relevant columns
    reducedModelMatrix <- modelMatrix[,idx,drop=FALSE]
    reducedBetaPriorVar <- betaPriorVar[idx]
    reducedLambda <- 1/reducedBetaPriorVar
    reducedModel <- fitNbinomGLMs(objectNZ, modelMatrix=reducedModelMatrix,
                                  lambda=reducedLambda,
                                  maxit=maxit, useOptim=useOptim, useQR=useQR)
  }
  
  attr(object,"betaPrior") <- betaPrior
  attr(object,"betaPriorVar") <- betaPriorVar
  attr(object,"modelMatrix") <- modelMatrix
  attr(object,"reducedModelMatrix") <- reducedModel$modelMatrix
  attr(object,"test") <- "LRT"

  # store mu in case the user did not call estimateDispersionsGeneEst
  dimnames(fullModel$mu) <- NULL
  assays(objectNZ)[["mu"]] <- fullModel$mu
  assays(object)[["mu"]] <- buildMatrixWithNARows(fullModel$mu, mcols(object)$allZero)

  H <- fullModel$hat_diagonals

  # calculate Cook's distance
  dispModelMatrix <- modelMatrix
  attr(object,"dispModelMatrix") <- dispModelMatrix
  cooks <- calculateCooksDistance(objectNZ, H, dispModelMatrix)
  
  # record maximum of Cook's
  maxCooks <- recordMaxCooks(design(object), colData(object), dispModelMatrix, cooks, nrow(objectNZ))

  # store Cook's distance for each sample
  assays(object)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)
  
  if (any(!fullModel$betaConv)) {
    if (!quiet) message(paste(sum(!fullModel$betaConv),"rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT"))
  }

  # calculate LRT statistic and p-values
  LRTStatistic <- (2 * (fullModel$logLike - reducedModel$logLike))
  LRTPvalue <- pchisq(LRTStatistic, df=df, lower.tail=FALSE)

  mleBetas <- if (betaPrior) {
    matrixToList(mleBetaMatrix)
  } else {
    NULL
  }
  
  # continue storing LRT results
  resultsList <- c(matrixToList(fullModel$betaMatrix),
                   matrixToList(fullModel$betaSE),
                   mleBetas,
                   list(LRTStatistic = LRTStatistic,
                        LRTPvalue = LRTPvalue,
                        fullBetaConv = fullModel$betaConv,
                        reducedBetaConv = reducedModel$betaConv,
                        betaIter = fullModel$betaIter,
                        deviance = -2 * fullModel$logLike,
                        maxCooks = maxCooks))
  LRTResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)

  modelComparison <- if (modelAsFormula) {
    paste0("'",paste(as.character(full),collapse=" "),
           "' vs '", paste(as.character(reduced),collapse=" "),"'")
  } else {
    "full vs reduced"
  }

  modelMatrixNames <- colnames(fullModel$betaMatrix)
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  lfcType <- if (attr(object,"betaPrior")) "MAP" else "MLE"
  coefInfo <- paste(paste0("log2 fold change (",lfcType,"):"),modelMatrixNamesSpaces)
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)
  mleInfo <- if (betaPrior) {
    gsub("_"," ",colnames(mleBetaMatrix))
  } else {
    NULL
  }
  statInfo <- paste("LRT statistic:",modelComparison)
  pvalInfo <- paste("LRT p-value:",modelComparison)

  mcols(LRTResults) <- DataFrame(type = rep("results",ncol(LRTResults)),
                                 description = c(coefInfo, seInfo, mleInfo,
                                   statInfo, pvalInfo, 
                                   "convergence of betas for full model",
                                   "convergence of betas for reduced model",
                                   "iterations for betas for full model",
                                   "deviance of the full model",
                                   "maximum Cook's distance for row"))
  mcols(object) <- cbind(mcols(object),LRTResults)
  return(object)
}


#' Replace outliers with trimmed mean
#'
#' Note that this function is called within \code{\link{DESeq}}, so is not
#' necessary to call on top of a \code{DESeq} call. See the \code{minReplicatesForReplace}
#' argument documented in \code{link{DESeq}}.
#' 
#' This function replaces outlier counts flagged by extreme Cook's distances,
#' as calculated by \code{\link{DESeq}}, \code{\link{nbinomWaldTest}}
#' or \code{\link{nbinomLRT}}, with values predicted by the trimmed mean
#' over all samples (and adjusted by size factor or normalization factor).
#' This function replaces the counts in the matrix returned by \code{counts(dds)}
#' and the Cook's distances in \code{assays(dds)[["cooks"]]}. Original counts are
#' preserved in \code{assays(dds)[["originalCounts"]]}.
#' 
#' The \code{\link{DESeq}} function calculates a diagnostic measure called
#' Cook's distance for every gene and every sample. The \code{\link{results}}
#' function then sets the p-values to \code{NA} for genes which contain
#' an outlying count as defined by a Cook's distance above a threshold.
#' With many degrees of freedom, i.e. many more samples than number of parameters to 
#' be estimated-- it might be undesirable to remove entire genes from the analysis
#' just because their data include a single count outlier.
#' An alternate strategy is to replace the outlier counts
#' with the trimmed mean over all samples, adjusted by the size factor or normalization
#' factor for that sample. The following simple function performs this replacement
#' for the user, for samples which have at least \code{minReplicates} number
#' of replicates (including that sample).
#' For more information on Cook's distance, please see the two
#' sections of the vignette: 'Dealing with count outliers' and 'Count outlier detection'.
#' 
#' @param object a DESeqDataSet object, which has already been processed by
#' either DESeq, nbinomWaldTest or nbinomLRT, and therefore contains a matrix
#' contained in \code{assays(dds)[["cooks"]]}. These are the Cook's distances which will
#' be used to define outlier counts.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from
#' each end of the normalized counts for a gene before the mean is computed
#' @param cooksCutoff the threshold for defining an outlier to be replaced.
#' Defaults to the .99 quantile of the F(p, m - p) distribution, where p is
#' the number of parameters and m is the number of samples.
#' @param minReplicates the minimum number of replicate samples necessary to consider
#' a sample eligible for replacement (including itself). Outlier counts will not be replaced
#' if the sample is in a cell which has less than minReplicates replicates.
#' @param whichSamples optional, a numeric or logical index to specify
#' which samples should have outliers replaced. if missing, this is determined using
#' minReplicates.
#'
#' @seealso \code{\link{DESeq}}
#'
#' @aliases replaceOutliersWithTrimmedMean
#' 
#' @return a DESeqDataSet with replaced counts in the slot returned by
#' \code{\link{counts}} and the original counts preserved in
#' \code{assays(dds)[["originalCounts"]]}
#' 
#' @export
replaceOutliers <- function(object, trim=.2, cooksCutoff, minReplicates=7, whichSamples) {
  if (is.null(attr(object,"modelMatrix")) | !("cooks" %in% assayNames(object))) {
    stop("first run DESeq, nbinomWaldTest, or nbinomLRT to identify outliers")
  }
  if (minReplicates < 3) {
    stop("at least 3 replicates are necessary in order to indentify a sample as a count outlier")
  }
  stopifnot(is.numeric(minReplicates) & length(minReplicates) == 1)
  p <- ncol(attr(object,"modelMatrix"))
  m <- ncol(object)
  if (m <= p) {
    assays(object)[["originalCounts"]] <- counts(object)
    return(object)
  }
  if (missing(cooksCutoff)) {
    cooksCutoff <- qf(.99, p, m - p)
  }
  idx <- which(assays(object)[["cooks"]] > cooksCutoff)
  mcols(object)$replace <- apply(assays(object)[["cooks"]], 1, function(row) any(row > cooksCutoff))
  mcols(mcols(object),use.names=TRUE)["replace",] <- DataFrame(type="intermediate",description="had counts replaced")
  trimBaseMean <- apply(counts(object,normalized=TRUE),1,mean,trim=trim)
  # build a matrix of counts based on the trimmed mean and the size factors
  replacementCounts <- if (!is.null(normalizationFactors(object))) {
    as.integer(matrix(rep(trimBaseMean,ncol(object)),ncol=ncol(object)) * 
               normalizationFactors(object))
  } else {
    as.integer(outer(trimBaseMean, sizeFactors(object), "*"))
  }
  
  # replace only those values which fall above the cutoff on Cook's distance
  newCounts <- counts(object)
  newCounts[idx] <- replacementCounts[idx]
  
  if (missing(whichSamples)) {
    whichSamples <- nOrMoreInCell(attr(object,"modelMatrix"), n = minReplicates)
  }
  stopifnot(is.logical(whichSamples))
  object$replaceable <- whichSamples
  mcols(colData(object),use.names=TRUE)["replaceable",] <- DataFrame(type="intermediate",
                         description="outliers can be replaced")
  assays(object)[["originalCounts"]] <- counts(object)
  if (sum(whichSamples) == 0) {
    return(object)
  }
  counts(object)[,whichSamples] <- newCounts[,whichSamples,drop=FALSE]
  object
}

#' @export
#' @rdname replaceOutliers
replaceOutliersWithTrimmedMean <- replaceOutliers


###########################################################
# unexported functons 
###########################################################


# Get base means and variances
#
# An internally used function to calculate the row means and variances
# from the normalized counts, which requires that \code{\link{estimateSizeFactors}}
# has already been called.  Adds these and a logical column if the row sums
# are zero to the mcols of the object.
#
# object a DESeqDataSet object
#
# return a DESeqDataSet object with columns baseMean
# and baseVar in the row metadata columns
getBaseMeansAndVariances <- function(object) {
  meanVarZero <- DataFrame(baseMean = unname(rowMeans(counts(object,normalized=TRUE))),
                           baseVar = unname(rowVars(counts(object,normalized=TRUE))),
                           allZero = unname(rowSums(counts(object)) == 0))
  mcols(meanVarZero) <- DataFrame(type = rep("intermediate",ncol(meanVarZero)),
                                  description = c("mean of normalized counts for all samples",
                                    "variance of normalized counts for all samples",
                                    "all counts for a gene are zero"))
  if (all(c("baseMean","baseVar","allZero") %in% names(mcols(object)))) {
      mcols(object)[c("baseMean","baseVar","allZero")] <- meanVarZero
  } else {
      mcols(object) <- cbind(mcols(object),meanVarZero)
  }
  return(object)
}

estimateNormFactors <- function(counts, normMatrix, locfunc=median, geoMeans, controlGenes) {
  sf <- estimateSizeFactorsForMatrix(counts / normMatrix, locfunc=locfunc, geoMeans=geoMeans, controlGenes=controlGenes)
  nf <- t( t(normMatrix) * sf )
  nf / exp(rowMeans(log(nf)))
}

# Estimate a parametric fit of dispersion to the mean intensity
parametricDispersionFit <- function( means, disps ) {
   coefs <- c( .1, 1 )
   iter <- 0
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which( (residuals > 1e-4) & (residuals < 15) )
      # check for glm convergence below to exit while-loop
      suppressWarnings({fit <- glm( disps[good] ~ I(1/means[good]),
         family=Gamma(link="identity"), start=coefs )})
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      if ( !all( coefs > 0 ) )
         stop(simpleError("parametric dispersion fit failed"))
      if ( ( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )  & fit$converged )
         break
      iter <- iter + 1
      if ( iter > 10 ) 
        stop(simpleError("dispersion fit did not converge"))
    }
   names( coefs ) <- c( "asymptDisp", "extraPois" )
   ans <- function(q) coefs[1] + coefs[2] / q
   attr( ans, "coefficients" ) <- coefs
   ans
}


# Local fit of dispersion to the mean intensity
# fitting is done on log dispersion, log mean scale
localDispersionFit <- function( means, disps, minDisp ) {
  if (all(disps < minDisp*10)) {
    return(rep(minDisp,length(disps)))
  }
  d <- data.frame(logDisps = log(disps), logMeans = log(means))
  fit <- locfit(logDisps ~ logMeans, data=d[disps >= minDisp*10,,drop=FALSE],
                weights = means[disps >= minDisp*10])
  dispFunction <- function(means) exp(predict(fit, data.frame(logMeans=log(means))))
  return(dispFunction)
}


# convenience function for testing the log likelihood
# for a count matrix, mu matrix and vector disp
nbinomLogLike <- function(counts, mu, disp) {
  rowSums(matrix(dnbinom(counts, mu=mu,size=1/disp,
                         log=TRUE),ncol=ncol(counts)))
}


# Unexported, low-level function for fitting negative binomial GLMs
#
# Users typically call \code{\link{nbinomWaldTest}} or \code{\link{nbinomLRT}}
# which calls this function to perform fitting.  These functions return
# a \code{\link{DESeqDataSet}} object with the appropriate columns
# added.  This function returns results as a list.
#
# object a DESeqDataSet
# modelMatrix the design matrix
# modelFormula a formula specifying how to construct the design matrix
# alpha_hat the dispersion parameter estimates
# lambda the 'ridge' term added for the penalized GLM on the log2 scale
# renameCols whether to give columns variable_B_vs_A style names
# betaTol control parameter: stop when the following is satisfied:
#   abs(dev - dev_old)/(abs(dev) + 0.1) < betaTol
# maxit control parameter: maximum number of iteration to allow for
#   convergence
# useOptim whether to use optim on rows which have not converged:
#   Fisher scoring is not ideal with multiple groups and sparse
#   count distributions
# useQR whether to use the QR decomposition on the design matrix X
# forceOptim whether to use optim on all rows
# warnNonposVar whether to warn about non positive variances,
#   for advanced users only running LRT without beta prior,
#   this might be desirable to be ignored.
#
# return a list of results, with coefficients and standard
# errors on the log2 scale
fitNbinomGLMs <- function(object, modelMatrix=NULL, modelFormula, alpha_hat, lambda,
                          renameCols=TRUE, betaTol=1e-8, maxit=100, useOptim=TRUE,
                          useQR=TRUE, forceOptim=FALSE, warnNonposVar=TRUE) {
  if (missing(modelFormula)) {
    modelFormula <- design(object)
  }
  if (is.null(modelMatrix)) {
    modelAsFormula <- TRUE
    modelMatrix <- model.matrix(modelFormula, data=colData(object))
  } else {
    modelAsFormula <- FALSE
  }

  stopifnot(all(colSums(abs(modelMatrix)) > 0))

  # rename columns, for use as columns in DataFrame
  # and to emphasize the reference level comparison
  modelMatrixNames <- colnames(modelMatrix)
  modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  modelMatrixNames <- make.names(modelMatrixNames)
  
  if (renameCols) {
    convertNames <- renameModelMatrixColumns(colData(object),
                                             modelFormula)
    convertNames <- convertNames[convertNames$from %in% modelMatrixNames,,drop=FALSE]
    modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
  }
  colnames(modelMatrix) <- modelMatrixNames
  
  normalizationFactors <- if (!is.null(normalizationFactors(object))) {
    normalizationFactors(object)
  } else { 
    matrix(rep(sizeFactors(object),each=nrow(object)),
           ncol=ncol(object))
  }
  
  if (missing(alpha_hat)) {
    alpha_hat <- dispersions(object)
  }

  if (length(alpha_hat) != nrow(object)) {
    stop("alpha_hat needs to be the same length as nrows(object)")
  }

  # set a wide prior for all coefficients
  if (missing(lambda)) {
    lambda <- rep(1e-6, ncol(modelMatrix))
  }
  
  # bypass the beta fitting if the model formula is only intercept and
  # the prior variance is large (1e6)
  # i.e., LRT with reduced ~ 1 and no beta prior
  justIntercept <- if (modelAsFormula) {
    modelFormula == formula(~ 1)
  } else {
    ncol(modelMatrix) == 1 & all(modelMatrix == 1)
  }
  if (justIntercept & all(lambda <= 1e-6)) {
      alpha <- alpha_hat
      betaConv <- rep(TRUE, nrow(object))
      betaIter <- rep(1,nrow(object))
      betaMatrix <- matrix(log2(mcols(object)$baseMean),ncol=1)
      mu <- normalizationFactors * as.numeric(2^betaMatrix)
      logLike <- rowSums(dnbinom(counts(object), mu=mu, size=1/alpha, log=TRUE))
      deviance <- -2 * logLike
      modelMatrix <- model.matrix(~ 1, colData(object))
      colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
      w <- (mu^-1 + alpha)^-1
      xtwx <- rowSums(w)
      sigma <- xtwx^-1
      betaSE <- matrix(log2(exp(1)) * sqrt(sigma),ncol=1)      
      hat_diagonals <- w * xtwx^-1;
      res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
                  betaSE = betaSE, mu = mu, betaIter = betaIter,
                  deviance = deviance,
                  modelMatrix=modelMatrix, 
                  nterms=1, hat_diagonals=hat_diagonals)
      return(res)
  }
  
  qrx <- qr(modelMatrix)
  # if full rank, estimate initial betas for IRLS below
  if (qrx$rank == ncol(modelMatrix)) {
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)
    y <- t(log(counts(object,normalized=TRUE) + .1))
    beta_mat <- t(solve(R, t(Q) %*% y))
  } else {
    if ("Intercept" %in% modelMatrixNames) {
      beta_mat <- matrix(0, ncol=ncol(modelMatrix), nrow=nrow(object))
      # use the natural log as fitBeta occurs in the natural log scale
      logBaseMean <- log(rowMeans(counts(object,normalized=TRUE)))
      beta_mat[,which(modelMatrixNames == "Intercept")] <- logBaseMean
    } else {
      beta_mat <- matrix(1, ncol=ncol(modelMatrix), nrow=nrow(object))
    }
  }
  
  # here we convert from the log2 scale of the betas
  # and the beta prior variance to the log scale
  # used in fitBeta.
  # so we divide by the square of the
  # conversion factor, log(2)
  lambdaLogScale <- lambda / log(2)^2

  betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix,
                            nfSEXP = normalizationFactors,
                            alpha_hatSEXP = alpha_hat,
                            beta_matSEXP = beta_mat,
                            lambdaSEXP = lambdaLogScale,
                            tolSEXP = betaTol, maxitSEXP = maxit,
                            useQRSEXP=useQR)
  mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
  dispersionVector <- rep(dispersions(object), times=ncol(object))
  logLike <- nbinomLogLike(counts(object), mu, dispersions(object))

  # test for stability
  rowStable <- apply(betaRes$beta_mat,1,function(row) sum(is.na(row))) == 0

  # test for positive variances
  rowVarPositive <- apply(betaRes$beta_var_mat,1,function(row) sum(row <= 0)) == 0
  
  # test for convergence, stability and positive variances
  betaConv <- betaRes$iter < maxit
  
  # here we transform the betaMatrix and betaSE to a log2 scale
  betaMatrix <- log2(exp(1))*betaRes$beta_mat
  colnames(betaMatrix) <- modelMatrixNames
  colnames(modelMatrix) <- modelMatrixNames
  # warn below regarding these rows with negative variance
  betaSE <- log2(exp(1))*sqrt(pmax(betaRes$beta_var_mat,0))
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)

  # switch based on whether we should also use optim
  # on rows which did not converge
  rowsForOptim <- if (useOptim) {
    which(!betaConv | !rowStable | !rowVarPositive)
  } else {
    which(!rowStable | !rowVarPositive)
  }
  
  if (forceOptim) {
    rowsForOptim <- seq_along(betaConv)
  }
  
  if (length(rowsForOptim) > 0) {
    # we use optim if didn't reach convergence with the IRLS code
    resOptim <- fitNbinomGLMsOptim(object,modelMatrix,lambda,
                                   rowsForOptim,rowStable,
                                   normalizationFactors,alpha_hat,
                                   betaMatrix,betaSE,betaConv,
                                   beta_mat,
                                   mu,logLike)
    betaMatrix <- resOptim$betaMatrix
    betaSE <- resOptim$betaSE
    betaConv <- resOptim$betaConv
    mu <- resOptim$mu
    logLike <- resOptim$logLike
  }

  stopifnot(!any(is.na(betaSE)))
  nNonposVar <- sum(rowSums(betaSE == 0) > 0)
  if (warnNonposVar & nNonposVar > 0) warning(nNonposVar,"rows had non-positive estimates of variance for coefficients")
  
  list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter,
       deviance = betaRes$deviance,
       modelMatrix=modelMatrix, 
       nterms=ncol(modelMatrix), hat_diagonals=betaRes$hat_diagonals)
}


# Fit dispersions for negative binomial GLM
#
# This function estimates the dispersion parameter (alpha) for negative binomial
# generalized linear models. The fitting is performed on the log scale.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# mu_hatSEXP n by m matrix, the expected mean values, given beta-hat
# log_alphaSEXP n length vector of initial guesses for log(alpha)
# log_alpha_prior_meanSEXP n length vector of the fitted values for log(alpha)
# log_alpha_prior_sigmasqSEXP a single numeric value for the variance of the prior
# min_log_alphaSEXP the minimum value of log alpha
# kappa_0SEXP a parameter used in calculting the initial proposal
#   for the backtracking search
#   initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
# use_priorSEXP boolean variable, whether to use a prior or just calculate the MLE
#
# return a list with elements: log_alpha, iter, iter_accept, last_change, initial_lp, intial_dlp, last_lp, last_dlp, last_d2lp
fitDispWrapper <- function (ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP,
                            log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP,
                            tolSEXP, maxitSEXP, use_priorSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDispWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitDisp, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  fitDisp(ySEXP=ySEXP, xSEXP=xSEXP, mu_hatSEXP=mu_hatSEXP,
          log_alphaSEXP=log_alphaSEXP, log_alpha_prior_meanSEXP=log_alpha_prior_meanSEXP,
          log_alpha_prior_sigmasqSEXP=log_alpha_prior_sigmasqSEXP,
          min_log_alphaSEXP=min_log_alphaSEXP, kappa_0SEXP=kappa_0SEXP,
          tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, use_priorSEXP=use_priorSEXP)
}


# Fit beta coefficients for negative binomial GLM
#
# This function estimates the coefficients (betas) for negative binomial generalized linear models.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# nfSEXP n by m matrix of normalization factors
# alpha_hatSEXP n length vector of the disperion estimates
# contrastSEXP a k length vector for a possible contrast
# beta_matSEXP n by k matrix of the initial estimates for the betas
# lambdaSEXP k length vector of the ridge values
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
# useQRSEXP whether to use QR decomposition
#
# Note: at this level the betas are on the natural log scale
fitBetaWrapper <- function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP,
                            beta_matSEXP, lambdaSEXP, tolSEXP, maxitSEXP, useQRSEXP) {
  if ( missing(contrastSEXP) ) {
    # contrast is not required, just give 1,0,0,...
    contrastSEXP <- c(1,rep(0,ncol(xSEXP)-1))
  }
  # test for any NAs in arguments
  arg.names <- names(formals(fitBetaWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitBeta, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  
  fitBeta(ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP, alpha_hatSEXP=alpha_hatSEXP,
          contrastSEXP=contrastSEXP, beta_matSEXP=beta_matSEXP,
          lambdaSEXP=lambdaSEXP, tolSEXP=tolSEXP, maxitSEXP=maxitSEXP,
          useQRSEXP=useQRSEXP)
}


# convenience function for building results tables
# out of a list and filling in NA rows
buildDataFrameWithNARows <- function(resultsList, NArows) {
  lengths <- sapply(resultsList,length)
  if (!all(lengths == lengths[1])) {
    stop("lengths of vectors in resultsList must be equal")
  }
  if (sum(!NArows) != lengths[1]) {
    stop("number of non-NA rows must be equal to lengths of vectors in resultsList")
  }
  if (sum(NArows) == 0) {
    return(DataFrame(resultsList))
  }
  dfFull <- DataFrame(lapply(resultsList, function(x) vector(mode(x), length(NArows))))
  dfFull[NArows,] <- NA
  dfFull[!NArows,] <- DataFrame(resultsList)
  dfFull
}

# convenience function for building larger matrices
# by filling in NA rows
buildMatrixWithNARows <- function(m, NARows) {
  mFull <- matrix(NA, ncol=ncol(m), nrow=length(NARows))
  mFull[!NARows,] <- m
  mFull
}

# convenience function for building larger matrices
# by filling in 0 rows
buildMatrixWithZeroRows <- function(m, zeroRows) {
  mFull <- matrix(0, ncol=ncol(m), nrow=length(zeroRows))
  mFull[!zeroRows,] <- m
  mFull
}

# convenience function for breaking up matrices
# by column and preserving column names
matrixToList <- function(m) {
  l <- split(m, col(m))
  names(l) <- colnames(m)
  l
}


# calculate a robust method of moments dispersion,
# in order to estimate the dispersion excluding
# individual outlier counts which would raise the variance estimate
robustMethodOfMomentsDisp <- function(object, modelMatrix) {
  cnts <- counts(object,normalized=TRUE)
  # if there are 3 or more replicates in any cell
  threeOrMore <- nOrMoreInCell(modelMatrix,n=3)
  v <- if (any(threeOrMore)) {
    cells <- apply(modelMatrix,1,paste0,collapse="")
    cells <- unname(factor(cells,levels=unique(cells)))
    levels(cells) <- seq_along(levels(cells))
    levelsThreeOrMore <- levels(cells)[table(cells) >= 3]
    idx <- cells %in% levelsThreeOrMore
    cntsSub <- cnts[,idx,drop=FALSE]
    cellsSub <- factor(cells[idx])
    trimmedCellVariance(cntsSub, cellsSub)
  } else {
    trimmedVariance(cnts)
  }
  m <- rowMeans(cnts)
  alpha <- ( v - m ) / m^2
  # cannot use the typical minDisp = 1e-8 here or else all counts in the same
  # group as the outlier count will get an extreme Cook's distance
  minDisp <- 0.04
  alpha <- pmax(alpha, minDisp)
  alpha
}

trimmedCellVariance <- function(cnts, cells) {
  # how much to trim at different n
  trimratio <- c(1/3, 1/4, 1/8)
  # returns an index for the vector above for three sample size bins
  trimfn <- function(n) as.integer(cut(n, breaks=c(0,3.5,23.5,Inf)))
  cellMeans <- matrix(sapply(levels(cells), function(lvl) {
    n <- sum(cells==lvl)
    apply(cnts[,cells==lvl,drop=FALSE],1,mean,trim=trimratio[trimfn(n)])
  }),
                      nrow=nrow(cnts))
  qmat <- cellMeans[,as.integer(cells),drop=FALSE]
  sqerror <- (cnts - qmat)^2
  varEst <- matrix(sapply(levels(cells), function(lvl) {
    n <- sum(cells==lvl)
    # scale due to trimming of large squares, by e.g. 1/mean(rnorm(1e6)^2,trim=1/8)
    scale.c <- c(2.04, 1.86, 1.51)[trimfn(n)]
    scale.c * apply(sqerror[,cells==lvl,drop=FALSE],1,mean,trim=trimratio[trimfn(n)])
  }),
                   nrow=nrow(sqerror))
  # take the max of variance estimates from cells
  # as one condition might have highly variable counts
  rowMax(varEst)
}

trimmedVariance <- function(x) {
  rm <-  apply(x,1,mean,trim=1/8)
  sqerror <- (x - rm)^2
  # scale due to trimming of large squares
  1.51 * apply(sqerror,1,mean,trim=1/8)
}

calculateCooksDistance <- function(object, H, modelMatrix) {
  p <- ncol(modelMatrix)
  dispersions <- robustMethodOfMomentsDisp(object, modelMatrix)
  V <- assays(object)[["mu"]] + dispersions * assays(object)[["mu"]]^2
  PearsonResSq <- (counts(object) - assays(object)[["mu"]])^2 / V
  cooks <- PearsonResSq / p  * H / (1 - H)^2
  cooks
}


# this function breaks out the logic for calculating the max Cook's distance:
# the samples over which max Cook's distance is calculated:
#
# Cook's distance is considered for those samples with 3 or more replicates per cell
#
# if m == p or there are no samples over which to calculate max Cook's, then give NA
recordMaxCooks <- function(design, colData, modelMatrix, cooks, numRow) {
    samplesForCooks <- nOrMoreInCell(modelMatrix, n=3)
    p <- ncol(modelMatrix)
    m <- nrow(modelMatrix)
    maxCooks <- if ((m > p) & any(samplesForCooks)) {
      apply(cooks[,samplesForCooks,drop=FALSE], 1, max)
    } else {
      rep(NA, numRow)
    }
    maxCooks
}

# for each sample in the model matrix,
# are there n or more replicates in the same cell
# (including that sample)
# so for a 2 x 3 comparison, the returned vector for n = 3 is:
# FALSE, FALSE, TRUE, TRUE, TRUE
nOrMoreInCell <- function(modelMatrix, n) {
  numEqual <- sapply(seq_len(nrow(modelMatrix)), function(i) {
    modelMatrixDiff <- t(t(modelMatrix) - modelMatrix[i,])
    sum(apply(modelMatrixDiff, 1, function(row) all(row == 0)))
  })
  numEqual >= n
}

# returns TRUE or FALSE if removing row would leave matrix full rank
## leaveOneOutFullRank <- function(modelMatrix) {
##   sapply(seq_len(nrow(modelMatrix)), function(i) qr(modelMatrix[-i,,drop=FALSE])$rank) == ncol(modelMatrix)
## }


# an unexported diagnostic function
# to retrieve the covariance matrix
# for the GLM coefficients of a single row
# only for standard model matrices
covarianceMatrix <- function(object, rowNumber) {
  if (attr(object, "modelMatrixType") != "standard")
    stop("only for standard model matrices")
  # convert coefficients to log scale
  coefColumns <- names(mcols(object))[grep("log2 fold change",mcols(mcols(object))$description)]
  beta <- log(2) * as.numeric(as.data.frame(mcols(object)[rowNumber,coefColumns,drop=FALSE]))
  x <- model.matrix(design(object), colData(object))
  y <- counts(object)[rowNumber,]
  sf <- sizeFactors(object)
  alpha <- dispersions(object)[rowNumber]
  mu.hat <- as.vector(sf * exp(x %*% beta))
  minmu <- 0.5
  mu.hat[mu.hat < minmu] <- minmu
  w <- diag(1/(1/mu.hat^2 * ( mu.hat + alpha * mu.hat^2 )))
  betaPriorVar <- attr(object,"betaPriorVar")
  ridge <- diag(1/(log(2)^2 * betaPriorVar))
  sigma <- solve(t(x) %*% w %*% x + ridge) %*% (t(x) %*% w %*% x) %*% t(solve(t(x) %*% w %*% x + ridge))
  # convert back to log2 scale
  sigmaLog2Scale <- log2(exp(1))^2 * sigma
  sigmaLog2Scale
}

getDesignFactors <- function(object) {
  design <- design(object)
  designVars <- all.vars(design)
  designVarsClass <- sapply(designVars, function(v) class(colData(object)[[v]]))
  designVars[designVarsClass == "factor"]
}

# looking at the values of x which are large
# in absolute value, find the zero-centered Normal distribution
# with the matching quantile, and return the variance
# of that Normal distribution
matchUpperQuantileForVariance <- function(x, upperQuantile=.05) {
  sdEst <- quantile(abs(x), 1 - upperQuantile) / qnorm(1 - upperQuantile/2)
  unname(sdEst)^2
}

matchWeightedUpperQuantileForVariance <- function(x, weights, upperQuantile=.05) {
  sdEst <- wtd.quantile(abs(x), weights=weights, 1 - upperQuantile, normwt=TRUE) / qnorm(1 - upperQuantile/2)
  unname(sdEst)^2
}





# this function calls fitNbinomGLMs() twice:
# 1 - without the beta prior, in order to calculate the
#     beta prior variance and hat matrix
# 2 - again but with the prior in order to get beta matrix and standard errors
fitGLMsWithPrior <- function(object, maxit, useOptim, useQR, betaPriorVar) {
  
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  modelMatrixType <- attr(object, "modelMatrixType")
  
  if (missing(betaPriorVar) | !("H" %in% assayNames(objectNZ))) {
    # first, fit the negative binomial GLM without a prior,
    # used to construct the prior variances
    # and for the hat matrix diagonals for calculating Cook's distance
    fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim, useQR=useQR,
                         renameCols = (modelMatrixType == "standard"))
    modelMatrix <- fit$modelMatrix
    modelMatrixNames <- colnames(modelMatrix)
    H <- fit$hat_diagonal
    betaMatrix <- fit$betaMatrix

    modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
    modelMatrixNames <- make.names(modelMatrixNames)
    colnames(betaMatrix) <- modelMatrixNames
    
    # save the MLE log fold changes for addMLE argument of results
    convertNames <- renameModelMatrixColumns(colData(object),
                                             design(objectNZ))
    convertNames <- convertNames[convertNames$from %in% modelMatrixNames,,drop=FALSE]
    modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
    mleBetaMatrix <- fit$betaMatrix
    colnames(mleBetaMatrix) <- paste0("MLE_",modelMatrixNames)

    # store for use in estimateBetaPriorVar below
    mcols(objectNZ) <- cbind(mcols(objectNZ), DataFrame(mleBetaMatrix))
  } else {
    # we can skip the first MLE fit because the
    # beta prior variance and hat matrix diagonals were provided
    modelMatrix <- model.matrix(design(objectNZ), colData(object))
    H <- assays(objectNZ)[["H"]]
    mleBetaMatrix <- as.matrix(mcols(objectNZ)[,grep("MLE_",names(mcols(objectNZ))),drop=FALSE])
  }
     
  if (missing(betaPriorVar)) {
    betaPriorVar <- estimateBetaPriorVar(objectNZ)
  } else {
    # else we are provided the prior variance:
    # check if the lambda is the correct length
    # given the design formula
    if (modelMatrixType == "expanded") {
      modelMatrix <- makeExpandedModelMatrix(objectNZ)
    }
    p <- ncol(modelMatrix)
    if (length(betaPriorVar) != p) {
      stop(paste("betaPriorVar should have length",p,"to match:",paste(colnames(modelMatrix),collapse=", ")))
    }
  }
  
  # refit the negative binomial GLM with a prior on betas
  if (any(betaPriorVar == 0)) {
    stop("beta prior variances are equal to zero for some variables")
  }
  lambda <- 1/betaPriorVar

  if (modelMatrixType == "standard") {
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda, maxit=maxit, useOptim=useOptim,
                         useQR=useQR)
    modelMatrix <- fit$modelMatrix
  } else {
    modelMatrix <- makeExpandedModelMatrix(objectNZ)
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda, maxit=maxit, useOptim=useOptim,
                         useQR=useQR, modelMatrix=modelMatrix, renameCols=FALSE)
  }

  res <- list(fit=fit, H=H, betaPriorVar=betaPriorVar,
              modelMatrix=modelMatrix, mleBetaMatrix=mleBetaMatrix)
  res
}


# breaking out the optim backup code from fitNbinomGLMs
fitNbinomGLMsOptim <- function(object,modelMatrix,lambda,
                               rowsForOptim,rowStable,
                               normalizationFactors,alpha_hat,
                               betaMatrix,betaSE,betaConv,
                               beta_mat,
                               mu,logLike) {
  scaleCols <- apply(modelMatrix,2,function(z) max(abs(z)))
  stopifnot(all(scaleCols > 0))
  x <- sweep(modelMatrix,2,scaleCols,"/")
  lambdaColScale <- lambda / scaleCols^2
  lambdaColScale <- ifelse(lambdaColScale == 0, 1e-6, lambdaColScale)
  lambdaLogScale <- lambda / log(2)^2
  lambdaLogScaleColScale <- lambdaLogScale / scaleCols^2
  large <- 30
  for (row in rowsForOptim) {
    betaRow <- if (rowStable[row] & all(abs(betaMatrix[row,]) < large)) {
      betaMatrix[row,] * scaleCols
    } else {
      beta_mat[row,] * scaleCols
    }
    nf <- normalizationFactors[row,]
    k <- counts(object)[row,]
    alpha <- alpha_hat[row]
    objectiveFn <- function(p) {
      mu_row <- as.numeric(nf * 2^(x %*% p))
      logLike <- sum(dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE))
      logPrior <- sum(dnorm(p,0,sqrt(1/lambdaColScale),log=TRUE))
      negLogPost <- -1 * (logLike + logPrior)
      if (is.finite(negLogPost)) negLogPost else 10^300
    }
    o <- optim(betaRow, objectiveFn, method="L-BFGS-B",lower=-large, upper=large)
    ridge <- if (length(lambdaLogScale) > 1) {
      diag(lambdaLogScaleColScale)
    } else {
      as.matrix(lambdaLogScaleColScale,ncol=1)
    }
    # if we converged, change betaConv to TRUE
    if (o$convergence == 0) {
      betaConv[row] <- TRUE
    }
    # with or without convergence, store the estimate from optim
    betaMatrix[row,] <- o$par / scaleCols
    # calculate the standard errors
    mu_row <- as.numeric(nf * 2^(x %*% o$par))
    minmu <- 0.5
    mu_row[mu_row < minmu] <- minmu
    w <- diag((mu_row^-1 + alpha)^-1)
    xtwx <- t(x) %*% w %*% x
    xtwxRidgeInv <- solve(xtwx + ridge)
    sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
    # warn below regarding these rows with negative variance
    betaSE[row,] <- log2(exp(1)) * sqrt(pmax(diag(sigma),0)) / scaleCols
    # store the new mu vector
    mu[row,] <- mu_row
    logLike[row] <- sum(dnbinom(k, mu=mu_row, size=1/alpha, log=TRUE))
  }
  return(list(betaMatrix=betaMatrix,betaSE=betaSE,
              betaConv=betaConv,
              mu=mu,logLike=logLike))
}


### DEPRECATED ###
# backup function in case dispersion doesn't converge
# this functionality is now implemented in Cpp below
fitDispInR <- function(y, x, mu, logAlphaPriorMean,
                       logAlphaPriorSigmaSq, usePrior) {
  disp <- numeric(nrow(y))
  # function to evaluate posterior
  logPost <- function(logAlpha) {
    alpha <- exp(logAlpha)
    w <- diag(1/(1/murow^2 * ( murow + alpha * murow^2 )))
    logLike <- sum(dnbinom(yrow, mu=murow, size=1/alpha, log=TRUE))
    coxReid <- -.5*(log(det(t(x) %*% w %*% x)))
    logPrior <- if (usePrior) {
      dnorm(logAlpha, logAlphaPriorMeanRow, sqrt(logAlphaPriorSigmaSq), log=TRUE)
    } else {
      0
    } 
    (logLike + coxReid + logPrior)
  }

  maxDisp <- max(10, ncol(y))
  s <- seq(from=log(1e-8),to=log(maxDisp),length=15)
  delta <- s[2] - s[1]
  # loop through rows
  for (i in seq_len(nrow(y))) {
    murow <- mu[i,]
    yrow <- y[i,]
    logAlphaPriorMeanRow <- logAlphaPriorMean[i]
    lpo <- sapply(s, logPost)
    sfine <- seq(from=s[which.max(lpo)]-delta, to=s[which.max(lpo)]+delta, length=15)
    lpofine <- sapply(sfine, logPost)
    disp[i] <- exp(sfine[which.max(lpofine)])
  }
  disp
}

# Fit dispersions by evaluating over grid
#
# This function estimates the dispersion parameter (alpha) for negative binomial
# generalized linear models. The fitting is performed on the log scale.
#
# ySEXP n by m matrix of counts
# xSEXP m by k design matrix
# mu_hatSEXP n by m matrix, the expected mean values, given beta-hat
# disp_gridSEXP the grid over which to estimate
# log_alpha_prior_meanSEXP n length vector of the fitted values for log(alpha)
# log_alpha_prior_sigmasqSEXP a single numeric value for the variance of the prior
# use_priorSEXP boolean variable, whether to use a prior or just calculate the MLE
#
# return a list with elements: 
fitDispGridWrapper <- function(y, x, mu, logAlphaPriorMean,
                               logAlphaPriorSigmaSq, usePrior) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDispGridWrapper))
  na.test <- sapply(mget(arg.names), function(x) any(is.na(x)))
  if (any(na.test)) stop(paste("in call to fitDispGridWrapper, the following arguments contain NA:",
                               paste(arg.names[na.test],collapse=", ")))
  minLogAlpha <- log(1e-8)
  maxLogAlpha <- log(max(10, ncol(y)))
  dispGrid <- seq(from=minLogAlpha, to=maxLogAlpha, length=15)
  logAlpha <- fitDispGrid(ySEXP=y, xSEXP=x, mu_hatSEXP=mu, disp_gridSEXP=dispGrid,
                          log_alpha_prior_meanSEXP=logAlphaPriorMean,
                          log_alpha_prior_sigmasqSEXP=logAlphaPriorSigmaSq,
                          use_priorSEXP=usePrior)$log_alpha
  exp(logAlpha)
}



# rough dispersion estimate using counts and fitted values
roughDispEstimate <- function(y, x) {

  # must be positive
  mu <- linearModelMu(y, x)
  mu <- matrix(pmax(1, mu), ncol=ncol(mu))
  
  m <- nrow(x)
  p <- ncol(x)

  # an alternate rough estimator with higher mean squared or absolute error
  # (rowSums( (y - mu)^2/(mu * (m - p)) ) - 1)/rowMeans(mu)
  
  # rough disp ests will be adjusted up to minDisp later
  est <- rowSums( ((y - mu)^2 - mu) / mu^2 ) / (m - p)
  pmax(est, 0)
}

momentsDispEstimate <- function(object) {
  xim <- if (!is.null(normalizationFactors(object))) {
    mean(1/colMeans(normalizationFactors(object)))
  } else {
    mean(1/sizeFactors(object))
  }
  bv <- mcols(object)$baseVar
  bm <- mcols(object)$baseMean
  (bv - xim*bm)/bm^2
}

# fast, rough estimation of means for rough dispersion estimation (above)
linearModelMu <- function(y, x) {
  qrx <- qr(x)
  Q <- qr.Q(qrx)
  Rinv <- solve(qr.R(qrx))
  hatmatrix <- x %*% Rinv %*% t(Q)
  t(hatmatrix %*% t(y))
}

# checks for LRT formulas, written as function to remove duplicate code
# in DESeq() and nbinomLRT()
checkLRT <- function(full, reduced) {
  reducedNotInFull <- !all.vars(reduced) %in% all.vars(full)
  if (any(reducedNotInFull)) {
    stop(paste("the following variables in the reduced formula not in the full formula:",
               paste(all.vars(reduced)[reducedNotInFull],collapse=", ")))
  }
}

# bulky code separated from DESeq()
refitWithoutOutliers <- function(object, test, betaPrior, full, reduced,
                                 quiet, minReplicatesForReplace, modelMatrix, modelMatrixType) {
  cooks <- assays(object)[["cooks"]]
  object <- replaceOutliers(object, minReplicates=minReplicatesForReplace)

  # refit without outliers, if there were any replacements
  nrefit <- sum(mcols(object)$replace, na.rm=TRUE)
  if ( nrefit > 0 ) {
    object <- getBaseMeansAndVariances(object)
    newAllZero <- which(mcols(object)$replace & mcols(object)$allZero)
  }
  # only refit if some of the replacements don't result in all zero counts
  # otherwise, these cases are handled by results()
  if ( nrefit > 0 && nrefit > length(newAllZero) ) {
    if (!quiet) message(paste("-- replacing outliers and refitting for", nrefit,"genes
-- DESeq argument 'minReplicatesForReplace' =",minReplicatesForReplace,"
-- original counts are preserved in counts(dds)"))
    
    # refit on those rows which had replacement
    refitReplace <- which(mcols(object)$replace & !mcols(object)$allZero)
    objectSub <- object[refitReplace,]
    intermediateOrResults <- which(mcols(mcols(objectSub))$type %in% c("intermediate","results"))
    mcols(objectSub) <- mcols(objectSub)[,-intermediateOrResults,drop=FALSE]

    # estimate gene-wise dispersion
    if (!quiet) message("estimating dispersions")
    objectSub <- estimateDispersionsGeneEst(objectSub, quiet=quiet, modelMatrix=modelMatrix)
    
    # need to redo fitted dispersion due to changes in base mean
    mcols(objectSub)$dispFit <- dispersionFunction(objectSub)(mcols(objectSub)$baseMean)
    mcols(mcols(objectSub),use.names=TRUE)["dispFit",] <- DataFrame(type="intermediate",
                             description="fitted values of dispersion")
    dispPriorVar <- attr( dispersionFunction(object), "dispPriorVar" )

    # estimate dispersion MAP
    objectSub <- estimateDispersionsMAP(objectSub, quiet=quiet,
                                        dispPriorVar=dispPriorVar, modelMatrix=modelMatrix)

    # fit GLM
    if (!quiet) message("fitting model and testing")
    if (test == "Wald") {
      betaPriorVar <- attr(object, "betaPriorVar")
      objectSub <- nbinomWaldTest(objectSub, betaPrior=betaPrior,
                                  betaPriorVar=betaPriorVar, quiet=quiet,
                                  modelMatrix=modelMatrix,
                                  modelMatrixType=modelMatrixType)
    } else if (test == "LRT") {
      if (!betaPrior) {
        objectSub <- nbinomLRT(objectSub, full=full, reduced=reduced, quiet=quiet)
      } else {
        betaPriorVar <- attr(object, "betaPriorVar")
        objectSub <- nbinomLRT(objectSub, full=full, reduced=reduced, betaPrior=betaPrior,
                               betaPriorVar=betaPriorVar, quiet=quiet)
      }
    }
    
    idx <- match(names(mcols(objectSub)), names(mcols(object)))
    mcols(object)[refitReplace, idx] <- mcols(objectSub)
    mcols(object)[newAllZero, mcols(mcols(object))$type == "results"] <- NA
    
    # continue to flag if some conditions have less than minReplicatesForReplace
    if (all(object$replaceable)) {
      mcols(object)$maxCooks <- NA
    } else {
      replaceCooks <- assays(object)[["cooks"]]
      replaceCooks[,object$replaceable] <- 0
      mcols(object)$maxCooks <- recordMaxCooks(design(object), colData(object),
                                               attr(object,"dispModelMatrix"), replaceCooks, nrow(object))
    }
  }
  
  if ( nrefit > 0 ) {
    # save the counts used for fitting as replaceCounts
    assays(object)[["replaceCounts"]] <- counts(object)
    assays(object)[["replaceCooks"]] <- assays(object)[["cooks"]]

    # preserve original counts and Cook's distances
    counts(object) <- assays(object)[["originalCounts"]]
    assays(object)[["cooks"]] <- cooks
    
    # no longer need this assay slot
    assays(object)[["originalCounts"]] <- NULL
  }
  
  object
}

sanitizeRowRanges <- function(object) {
  if (is.null(mcols(mcols(object)))) {
    mcols(mcols(object)) <- DataFrame(type=rep("input",ncol(mcols(object))),
                                      description=character(ncol(mcols(object))))
  }
  class(mcols(mcols(object))$type) <- "character"
  class(mcols(mcols(object))$description) <- "character"
  mcols(mcols(object))$type[ is.na(mcols(mcols(object))$type) ] <- ""
  mcols(mcols(object))$description[ is.na(mcols(mcols(object))$description) ] <- ""
  object
}

sanitizeColData <- function(object) {
  if (is.null(mcols(colData(object)))) {
    mcols(colData(object)) <- DataFrame(type=rep("input",ncol(colData(object))),
                                        description=character(ncol(colData(object))))
  }
  class(mcols(colData(object))$type) <- "character"
  class(mcols(colData(object))$description) <- "character"
  mcols(colData(object))$type[ is.na(mcols(colData(object))$type) ] <- ""
  mcols(colData(object))$description[ is.na(mcols(colData(object))$description) ] <- ""
  object
}

estimateSizeFactorsIterate <- function(object, niter=10, Q=0.05) {
  design(object) <- ~ 1
  sf <- rep(1, ncol(object))
  idx <- rowSums(counts(object)) > 0
  cts <- counts(object)[idx,]
  for (i in seq_len(niter)) {
    sizeFactors(object) <- sf
    object <- estimateDispersions(object, fitType="mean", quiet=TRUE)
    q <- t(t(assays(object)[["mu"]])/sf)[idx,]
    disps <- dispersions(object)[idx]
    sf.old <- sf
    fn <- function(p) {
      sf <- exp(p - mean(p))
      mu.new <- t(t(q) * sf)
      ll <- matrix(dnbinom(cts, mu=mu.new, size=1/disps, log=TRUE), ncol=ncol(cts))
      gene.ll <- rowSums(ll)
      sum(gene.ll[ gene.ll > quantile(gene.ll, Q) ])
    }
    res <- optim(log(sf.old), fn, control=list(fnscale=-1), method="L-BFGS-B")
    if (res$convergence != 0) {
      stop("iterative size factor normalization did not converge within an iteration")
    }
    sf <- exp(res$par - mean(res$par))
    # loop more than once, and test for convergence
    if (i > 1 & sum((log(sf.old) - log(sf))^2) < 1e-4) {
      break
    } else {
      if (i == niter) {
        stop("iterative size factor normalization did not converge")
      }
    }
  }
  sf
}

checkFullRank <- function(modelMatrix) {
  if (qr(modelMatrix)$rank < ncol(modelMatrix)) {
    if (any(apply(modelMatrix, 2, function(col) all(col == 0)))) {
      stop("the model matrix is not full rank, so the model cannot be fit as specified.
  Levels or combinations of levels without any samples have resulted in
  column(s) of zeros in the model matrix.

  Please read the vignette section 'Model matrix not full rank':

  vignette('DESeq2')")
    } else {
      stop("the model matrix is not full rank, so the model cannot be fit as specified.
  One or more variables or interaction terms in the design formula are linear
  combinations of the others and must be removed.

  Please read the vignette section 'Model matrix not full rank':

  vignette('DESeq2')")
    }
  }
}

designAndArgChecker <- function(object, betaPrior) {
  termsOrder <- attr(terms.formula(design(object)),"order")
  hasIntercept <- attr(terms(design(object)),"intercept") == 1
  interactionPresent <- any(termsOrder > 1)
  if (betaPrior & !hasIntercept) {
    stop("betaPrior=TRUE can only be used if the design has an intercept.
  if specifying + 0 in the design formula, use betaPrior=FALSE")
  }
  if (betaPrior & interactionPresent) {
    stop("betaPrior=FALSE should be used for designs with interactions")
  }

  design <- design(object)
  designVars <- all.vars(design)
  if (length(designVars) > 0) {
    if (any(sapply(designVars, function(v) any(is.na(colData(object)[[v]]))))) {
      stop("variables in the design formula cannot have NA values")
    }
    designFactors <- designVars[sapply(designVars, function(v) is(colData(object)[[v]], "factor"))]
    if (length(designFactors) > 0 && any(sapply(designFactors,function(v) any(table(colData(object)[[v]]) == 0)))) {
      stop("factors in design formula must have samples for each level.
  this error can arise when subsetting a DESeqDataSet, in which
  all the samples for one or more levels of a factor in the design were removed.
  if this was intentional, use droplevels() to remove these levels, e.g.:

  dds$condition <- droplevels(dds$condition)
")
    }
    if (any(sapply(designVars, function(v) is(colData(object)[[v]], "ordered")))) {
      stop("the design contains an ordered factor. The internal steps
that estimate the beta prior variance and produce resultsNames
do not work on ordered factors. You should instead use model.matrix()
and then provide your custom matrix to 'full' argument of DESeq().
(You should also provide a matrix to 'reduced' for test='LRT'.)")
    }
  }
}

