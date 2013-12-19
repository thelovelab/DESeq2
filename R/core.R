#' Differential expression analysis based on the negative binomial distribution
#'
#' This function performs a default analysis through the steps:
#' \enumerate{
#' \item estimation of size factors: \code{\link{estimateSizeFactors}}
#' \item estimation of dispersion: \code{\link{estimateDispersions}}
#' \item negative binomial GLM fitting and Wald statistics: \code{\link{nbinomWaldTest}}
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
#' \deqn{ \mu_{ij} = s_j q_{ij} }{ mu_ij = s_j * q_ij }
#' \deqn{ \log_2(q_{ij}) = x_{j.} \beta_i }{ log2(q_ij) = x_j. * beta_i }
#'
#' where counts \eqn{K_{ij}}{K_ij} for gene i, sample j are modeled using
#' a negative binomial distribution with fitted mean \eqn{\mu_{ij}}{mu_ij}
#' and a gene-specific dispersion parameter \eqn{\alpha_i}{alpha_i}.
#' The fitted mean is composed of a sample-specific size factor
#' \eqn{s_j}{s_j} and a parameter \eqn{q_{ij}}{q_ij} proportional to the
#' expected true concentration of fragments for sample j.
#' The coefficients \eqn{\beta_i}{beta_i} give the log2 fold changes for gene i for each
#' column of the model matrix \eqn{X}{X}.
#' The sample-specific size factors can be replaced by
#' gene-specific normalization factors for each sample using
#' \code{\link{normalizationFactors}}.  For details on the fitting of the log2
#' fold changes and calculation of p-values see \code{\link{nbinomWaldTest}}
#' (or \code{\link{nbinomLRT}} if using \code{test="LRT"}).
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
#' This replacement is performed by the \code{\link{replaceOutliersWithTrimmedMean}}
#' function. This default behavior helps to prevent filtering genes
#' based on Cook's distance when there are many degrees of freedom.
#' See \code{\link{results}} for more information about filtering using
#' Cook's distance, and the sections of the vignette. Original counts are
#' kept in the slot returned by \code{\link{counts}}, while replacement
#' counts which were used for testing are kept in
#' \code{assays(object)[["replaceCounts"]]}.
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
#' the non-intercept coefficients (Tikhonov/ridge regularization)
#' See \code{\link{nbinomWaldTest}} for description. Only used
#' for the Wald test.
#' @param full the full model formula, this should be the formula in
#' \code{design(object)}, only used by the likelihood ratio test
#' @param reduced a reduced formula to compare against, e.g.
#' the full model with a term or terms of interest removed,
#' only used by the likelihood ratio test
#' @param quiet whether to print messages at each step
#' @param minReplicatesForReplace the minimum number of replicates required
#' in order to use \code{\link{replaceOutliersWithTrimmedMean}} on a
#' sample. If there are samples with so many replicates, the model will
#' be refit after these replacing outliers, flagged by Cook's distance.
#' Set to \code{Inf} in order to never replace outliers.
#' 
#' @author Michael Love
#'
#' @references Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 11 (2010) R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#'
#' @import BiocGenerics GenomicRanges IRanges Rcpp RcppArmadillo methods
#' @importFrom locfit locfit
#' @importFrom genefilter rowVars filtered_p
#' 
#' @useDynLib DESeq2
#'
#' @seealso \code{\link{nbinomWaldTest}}, \code{\link{nbinomLRT}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(betaSD=1)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
#' resLRT <- results(ddsLRT)
#'
#' @export
DESeq <- function(object, test=c("Wald","LRT"),
                  fitType=c("parametric","local","mean"), betaPrior,
                  full=design(object), reduced, quiet=FALSE,
                  minReplicatesForReplace=7) {
  if (missing(test)) {
    test <- match.arg(test, choices=c("Wald","LRT"))
  }
  if (missing(betaPrior)) {
    betaPrior <- test == "Wald"
  }
  attr(object, "betaPrior") <- betaPrior
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
  if (!quiet) message("estimating dispersions")
  object <- estimateDispersions(object, fitType=fitType, quiet=quiet)
  if (!quiet) message("fitting model and testing")
  if (test == "Wald") {
    object <- nbinomWaldTest(object, betaPrior=betaPrior, quiet=quiet)                             
  } else if (test == "LRT") {
    object <- nbinomLRT(object, full=full, reduced=reduced, quiet=quiet)
  }
  # refit without outliers
  if (any(nOrMoreInCell(attr(object,"modelMatrix"),minReplicatesForReplace))) {
    if (!quiet) message(paste("-- refitting without outliers: samples with >=",
                              minReplicatesForReplace,"replicates
-- original counts retained in assays(dds,'originalCounts')"))
    object <- replaceOutliersWithTrimmedMean(object,minReplicates=minReplicatesForReplace)
    if (!quiet) message("estimating dispersions")
    object <- estimateDispersions(object, fitType=fitType, quiet=quiet)
    if (!quiet) message("fitting model and testing")
    if (test == "Wald") {
      object <- nbinomWaldTest(object, betaPrior=betaPrior, quiet=quiet)
    } else if (test == "LRT") {
      object <- nbinomLRT(object, full=full, reduced=reduced, quiet=quiet)
    }
    assays(object)[["replaceCounts"]] <- counts(object)
    counts(object) <- assays(object)[["originalCounts"]]
  }
  object
}

#' Make a simulated DESeqDataSet
#'
#' Constructs a simulated dataset of negative binomial data from
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
  colData <- DataFrame(sample = paste("sample",1:m,sep=""),
                       condition=factor(rep(c("A","B"),times=c(ceiling(m/2),floor(m/2)))))
  if (m > 1) {
    x <- model.matrix(~ colData$condition)
  } else {
    x <- cbind(rep(1,m),rep(0,m))
  }
  mu <- 2^(t(x %*% t(beta))) * rep(sizeFactors, each=n)
  countData <- matrix(rnbinom(m*n, mu=mu, size=1/dispersion), ncol=m)
  rownames(colData) <- colData$sample
  rowData <- GRanges("1",IRanges(start=(1:n - 1) * 100 + 1,width=100))
  names(rowData) <- paste0("feature",1:n)
  if (m > 1) {
    designFormula <- formula(~ condition)
  } else {
    designFormula <- formula(~ 1)
  }
  object <- DESeqDataSetFromMatrix(countData = countData,
                                   colData = colData,
                                   design = designFormula,
                                   rowData = rowData)
  trueVals <- DataFrame(trueIntercept = beta[,1],
                        trueBeta = beta[,2],
                        trueDisp = dispersion)
  mcols(trueVals) <- DataFrame(type=rep("input",ncol(trueVals)),
                               description=c("simulated intercept values",
                                 "simulated beta values",
                                 "simulated dispersion values"))
  mcols(object) <- cbind(mcols(object),trueVals)
  # cleaning up environment from formula above
  objNames <- ls()
  objNames <- objNames[objNames != "object"]
  rm(list=objNames)
  rm("objNames")
  object
}


#' Low-level function to estimate size factors with robust regression.
#' 
#' Given a matrix or data frame of count data, this function estimates the size
#' factors as follows: Each column is divided by the geometric means of the
#' rows. The median (or, if requested, another location estimator) of these
#' ratios (skipping the genes with a geometric mean of zero) is used as the size
#' factor for this column. Typically, you will not call this function directly, but use
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
estimateSizeFactorsForMatrix <- function( counts, locfunc = median, geoMeans)
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
  apply( counts, 2, function(cnts)
        exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & (cnts > 0)])))
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
#' for calculations in log scale, one decade above this value is used
#' as a test for inclusion in mean-dispersion fitting
#' @param kappa_0 control parameter used in setting the initial proposal
#' in backtracking search, higher kappa_0 results in larger steps
#' @param dispTol control parameter to test for convergence of log dispersion,
#' stop when increase in log posterior is less than dispTol
#' @param maxit control parameter: maximum number of iterations to allow for convergence
#' @param quiet whether to print messages at each step
#' @param modelMatrixType either "standard" or "expanded" for which
#' model matrix will be used later by \code{\link{nbinomWaldTest}}.
#' in case of "expanded" a full-rank model matrix is provided inside
#' this function which will not change from releveling.
#'
#' @return a DESeqDataSet with gene-wise, fitted, or final MAP
#' dispersion estimates in the metadata columns of the object.
#'
#' @aliases estimateDispersionsGeneEst estimateDispersionsFit estimateDispersionsMAP
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
#' @seealso \code{\link{estimateDispersions}}
#'
#' @export
estimateDispersionsGeneEst <- function(object, minDisp=1e-8, kappa_0=1,
                                       dispTol=1e-6, maxit=100, quiet=FALSE,
                                       modelMatrixType) {
  if ("dispGeneEst" %in% names(mcols(object))) {
    if (!quiet) message("you had estimated gene-wise dispersions, removing these")
    mcols(object) <- mcols(object)[,!names(mcols(object))  == "dispGeneEst"]
  }
  stopifnot(length(minDisp) == 1)
  stopifnot(length(kappa_0) == 1)
  stopifnot(length(dispTol) == 1)
  stopifnot(length(maxit) == 1)
  if (log(minDisp/10) <= -30) {
    stop("for computational stability, log(minDisp/10) should be above -30")
  }
  if (missing(modelMatrixType)) {
    betaPriorOrNotSpecified <- is.null(attr(object,"betaPrior")) || attr(object,"betaPrior")
    if (factorPresentThreeOrMoreLevels(object) & betaPriorOrNotSpecified) {
      modelMatrixType <- "expanded"
    } else {
      modelMatrixType <- "standard"
    }
  }
  object <- getBaseMeansAndVariances(object)
  if (!is.null(normalizationFactors(object))) {
    xim <- mean(1/colMeans(normalizationFactors(object)))
  } else {
    xim <- mean(1/sizeFactors(object))
  }
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]
  bv <- mcols(objectNZ)$baseVar
  bm <- mcols(objectNZ)$baseMean
  # this rough dispersion estimate (alpha_hat)
  # is for estimating beta and mu
  # and for the initial starting point for line search
  alpha_hat <- pmax(minDisp, (bv - xim*bm)/bm^2)

  if (modelMatrixType == "standard") {
    # fitNbinomGLMs returns mu and modelMatrix
    fit <- fitNbinomGLMs(objectNZ, alpha_hat=alpha_hat)
  } else {
    mm <- makeReleveledModelMatrix(object)
    fit <- fitNbinomGLMs(objectNZ, alpha_hat=alpha_hat, modelMatrix=mm)
  }
  
  # use of kappa_0 in backtracking search
  # initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
  # use log(minDisp/10) to stop if dispersions going to -infinity
  dispRes <- fitDisp(ySEXP = counts(objectNZ), xSEXP = fit$modelMatrix, mu_hatSEXP = fit$mu,
                     log_alphaSEXP = log(alpha_hat), log_alpha_prior_meanSEXP = log(alpha_hat),
                     log_alpha_prior_sigmasqSEXP = 1, min_log_alphaSEXP = log(minDisp/10),
                     kappa_0SEXP = kappa_0, tolSEXP = dispTol,
                     maxitSEXP = maxit, use_priorSEXP = FALSE)

  if (mean(dispRes$iter < maxit) < .5) {
    warning("in calling estimateDispersionsGeneEst, less than 50% of gene-wise estimates converged. Use larger maxit argument with estimateDispersions")
  }
  
  # dont accept moves if the log posterior did not
  # increase by more than one millionth,
  # and set the small estimates to the minimum dispersion
  dispGeneEst <- exp(dispRes$log_alpha)
  noIncrease <- dispRes$last_lp < dispRes$initial_lp + abs(dispRes$initial_lp)/1e6
  dispGeneEst[which(noIncrease)] <- alpha_hat[which(noIncrease)]
  dispGeneEst <- pmax(dispGeneEst, minDisp)
  dispGeneEstConv <- dispRes$iter < maxit
  
  dispDataFrame <- buildDataFrameWithNARows(list(dispGeneEst=dispGeneEst,
                                                 dispGeneEstConv=dispGeneEstConv),
                                            mcols(object)$allZero)
  mcols(dispDataFrame) <- DataFrame(type=rep("intermediate",ncol(dispDataFrame)),
                                    description=c("gene-wise estimates of dispersion",
                                      "gene-wise dispersion estimate convergence"))
  mcols(object) <- cbind(mcols(object), dispDataFrame)

  assays(object)[["mu"]] <- buildMatrixWithNARows(fit$mu, mcols(object)$allZero)
  
  return(object)
}

#' @rdname estimateDispersionsGeneEst
#' @export
estimateDispersionsFit <- function(object,fitType=c("parametric","local","mean"),
                                   minDisp=1e-8, quiet=FALSE) {
  if ("dispFit" %in% names(mcols(object))) {
    if (!quiet) message("you had estimated fitted dispersions, removing these")
    mcols(object) <- mcols(object)[,!names(mcols(object)) == "dispFit"]
  }
  objectNZ <- object[!mcols(object)$allZero,]
  useForFit <- mcols(objectNZ)$dispGeneEstConv

  # take the first fitType
  fitType <- fitType[1]
  stopifnot(length(fitType)==1)
  stopifnot(length(minDisp)==1)
  if (fitType == "parametric") {
    trial <- try(dispFunction <- parametricDispersionFit(mcols(objectNZ)$baseMean[useForFit],
                                                         mcols(objectNZ)$dispGeneEst[useForFit]),
                 silent=TRUE)
    if (!inherits(trial,"try-error")) {
      dispFit <- dispFunction(mcols(objectNZ)$baseMean)
    } else {
      warning("the parametric fit of dispersion estimates over the mean of counts
failed, which occurs when the trend is not well captured by the
function y = a/x + b. A local regression fit is automatically performed,
and the analysis can continue. You can specify fitType='local' or 'mean'
to avoid this message if re-running the same data.
When using local regression fit, the user should examine plotDispEsts(dds)
to make sure the fitted line is not sharply curving up or down based on
the position of individual points.")
      fitType <- "local"
    }
  }
  if (fitType == "local") {
    dispFunction <- localDispersionFit(means = mcols(objectNZ)$baseMean[useForFit],
                                       disps = mcols(objectNZ)$dispGeneEst[useForFit],
                                       minDisp = minDisp)
    dispFit <- dispFunction(mcols(objectNZ)$baseMean)
  }
  if (fitType == "mean") {
    useForMean <- mcols(objectNZ)$dispGeneEst > 10*minDisp
    meanDisp <- mean(mcols(objectNZ)$dispGeneEst[useForMean],na.rm=TRUE,trim=.05)
    dispFunction <- function(means) meanDisp
    dispFit <- rep(meanDisp,nrow(objectNZ))
  }
  if (!(fitType %in% c("parametric","local","mean"))) {
    stop("unknown fitType")
  }

  # store the dispersion function
  attr( dispFunction, "fitType" ) <- fitType
  dispersionFunction(object) <- dispFunction
  
  dispDataFrame <- buildDataFrameWithNARows(list(dispFit=dispFit),
                                            mcols(object)$allZero)
  mcols(dispDataFrame) <- DataFrame(type="intermediate",
                                    description="fitted values of dispersion")
  mcols(object) <- cbind(mcols(object), dispDataFrame)
  return(object)
}

#' @rdname estimateDispersionsGeneEst
#' @export
estimateDispersionsMAP <- function(object, outlierSD=2, dispPriorVar,
                                   minDisp=1e-8, kappa_0=1, dispTol=1e-6,
                                   maxit=100, quiet=FALSE) {
  stopifnot(length(outlierSD)==1)
  stopifnot(length(minDisp)==1)
  stopifnot(length(kappa_0)==1)
  stopifnot(length(dispTol)==1)
  stopifnot(length(maxit)==1)
  if ("dispersion" %in% names(mcols(object))) {
    if (!quiet) message("you had estimated dispersions, removing these")
    mcols(object) <- mcols(object)[,!names(mcols(object))  %in% c("dispersion","dispIter","dispIterAccept","dispConv")]
  }
 
  modelMatrix <- model.matrix(design(object), data=as.data.frame(colData(object)))  
  objectNZ <- object[!mcols(object)$allZero,]

  useNotMinDisp <- mcols(objectNZ)$dispGeneEst >= minDisp*10
  if (sum(useNotMinDisp,na.rm=TRUE) == 0) {
    warning(paste0("all genes have dispersion estimates < ",minDisp*10,
                   ", returning disp = ",minDisp*10))
    resultsList <- list(dispersion = minDisp*10)
    dispDataFrame <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
    mcols(dispDataFrame) <- DataFrame(type="intermediate",
                                      description="final estimates of dispersion")
    mcols(object) <- cbind(mcols(object), dispDataFrame)
    return(object)
  }
  
  # estimate the variance of the distribution of the
  # log dispersion estimates around the fitted value
  dispResiduals <- log(mcols(objectNZ)$dispGeneEst) - log(mcols(objectNZ)$dispFit)

  useForPrior <- useNotMinDisp
  if (sum(useForPrior,na.rm=TRUE) == 0) {
    stop("no data found which is greater than minDisp, within quants, and converged in gene-wise estimates")
  }

  varLogDispEsts <- varLogDispEstsAll <- mad(dispResiduals[useForPrior],na.rm=TRUE)^2
  
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
    obsDist <- dispResiduals[useForPrior]
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
    varLogDispEsts <- argminKL + expVarLogDisp
    # finally, restore the .Random.seed if it existed beforehand
    if (exists("oldRandomSeed")) {
      .Random.seed <<- oldRandomSeed
    }
  }

  attr( dispersionFunction(object), "varLogDispEsts" ) <- varLogDispEsts

  # estimate the expected sampling variance of the log estimates
  # Var(log(cX)) = Var(log(X))
  # X ~ chi-squared with m - p degrees of freedom
  if (m > p) {
    expVarLogDisp <- trigamma((m - p)/2)
    attr( dispersionFunction(object), "expVarLogDisp" ) <- expVarLogDisp
    # set the variance of the prior using these two estimates
    # with a minimum of .25
    dispPriorVarCalc <- pmax((varLogDispEsts - expVarLogDisp), .25)
  } else {
    # we have m = p, so do not try to substract sampling variance
    dispPriorVarCalc <- varLogDispEsts
  }
  

  # fill in the calculated dispersion prior variance
  if (missing(dispPriorVar)) {
    dispPriorVar <- dispPriorVarCalc
  }
  stopifnot(length(dispPriorVar)==1)
  attr( dispersionFunction(object), "dispPriorVar" ) <- dispPriorVar
  
  # set prior variance for fitting dispersion
  log_alpha_prior_sigmasq <- dispPriorVar

  # get previously calculated mu
  mu <- assays(objectNZ)[["mu"]]

  # start fitting at gene estimate unless the points are one decade
  # below the fitted line, then start at fitted line
  dispInit <- ifelse(mcols(objectNZ)$dispGeneEst >  0.1 * mcols(objectNZ)$dispFit,
                     mcols(objectNZ)$dispGeneEst,
                     mcols(objectNZ)$dispFit)

  # if any missing values, fill in the fitted value to initialize
  dispInit[is.na(dispInit)] <- mcols(objectNZ)$dispFit[is.na(dispInit)]
  
  # run with prior
  dispResMAP <- fitDisp(ySEXP = counts(objectNZ), xSEXP = modelMatrix, mu_hatSEXP = mu,
                        log_alphaSEXP = log(dispInit),
                        log_alpha_prior_meanSEXP = log(mcols(objectNZ)$dispFit),
                        log_alpha_prior_sigmasqSEXP = log_alpha_prior_sigmasq,
                        min_log_alphaSEXP = log(minDisp/10),
                        kappa_0SEXP = kappa_0, tolSEXP = dispTol,
                        maxitSEXP = maxit, use_priorSEXP = TRUE)

  # prepare dispersions for storage in mcols(object)
  dispersionFinal <- dispMAP <- exp(dispResMAP$log_alpha) 
    
  # detect outliers which have gene-wise estimates
  # outlierSD * standard deviation of log gene-wise estimates
  # above the fitted mean (prior mean)
  # and keep the original gene-est value for these.
  # Note: we use the variance of log dispersions estimates
  # from all the genes, not only those from below
  dispOutlier <- log(mcols(objectNZ)$dispGeneEst) >
                 log(mcols(objectNZ)$dispFit) +
                 outlierSD * sqrt(varLogDispEstsAll)
  dispOutlier[is.na(dispOutlier)] <- FALSE
  dispersionFinal[dispOutlier] <- mcols(objectNZ)$dispGeneEst[dispOutlier]
 
  resultsList <- list(dispersion = dispersionFinal,
                      dispIter = dispResMAP$iter,
                      dispConv = (dispResMAP$iter < maxit),
                      dispOutlier = dispOutlier,
                      dispMAP = dispMAP)

  numnonconv = sum(!resultsList$dispConv)
  if ((numnonconv>0) && !quiet) 
    message(sprintf("%d row%s did not converge in dispersion, labelled in 'mcols(object)$dispConv'. Try using a larger value for the 'maxit' argument of 'estimateDispersions'.\n", numnonconv, if(numnonconv>1) "s" else ""))
  
  dispDataFrame <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
  mcols(dispDataFrame) <- DataFrame(type=rep("intermediate",ncol(dispDataFrame)),
                                    description=c("final estimate of dispersion",
                                      "number of iterations",
                                      "convergence of final estimate",
                                      "dispersion flagged as outlier",
                                      "maximum a posteriori estimate"))

  mcols(object) <- cbind(mcols(object), dispDataFrame)
  return(object)
}



#' Wald test for the GLM coefficients
#' 
#' This function tests for significance of coefficients in a negative
#' binomial GLM, using previously calculated \code{\link{sizeFactors}}
#' (or \code{\link{normalizationFactors}})
#' and dispersion estimates.  See \code{\link{DESeq}} for the GLM formula.
#' 
#' The fitting proceeds as follows: standard maximum likelihood estimates
#' for GLM coefficients are calculated; a zero-mean normal prior distribution
#' is assumed; the variance of the prior distribution for each
#' non-intercept coefficient is calculated as the mean squared
#' maximum likelihood estimates over the genes which do not
#' contain zeros for some condition;
#' the final coefficients are then maximum a posteriori estimates
#' (using Tikhonov/ridge regularization) using this prior.
#' The use of a prior has little effect on genes with high counts and helps to
#' moderate the large spread in coefficients for genes with low counts.
#'
#' For calculating Wald test p-values, the coefficients are scaled by their
#' standard errors and then compared to a normal distribution. From
#' examination of Wald statistics for real datasets, the effect of the
#' prior on dispersion estimates results in a Wald statistic
#' distribution which is approximately normal.
#'
#' When interaction terms are present, the prior on log fold changes
#' the calculated beta prior variance will only be used for the interaction
#' terms (non-interaction terms receive a wide prior variance of 1e6).
#' In the case of interaction terms and factors with 3 or more
#' levels present in the design formula, a moderately wide prior
#' variance of 1e3 will be used on non-interaction terms, to allow for
#' convergence to the maximum of the posterior.
#' 
#' The Wald test can be replaced with the \code{\link{nbinomLRT}}
#' for an alternative test of significance.
#'
#' @param object a DESeqDataSet
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients (Tikhonov/ridge regularization)
#' @param betaPriorVar a vector with length equal to the number of
#' model terms including the intercept.
#  betaPriorVar gives the variance of the prior on the sample betas,
#' which if missing is estimated from the rows which do not have any
#' zeros
#' @param modelMatrixType either "standard" or "expanded", which describe
#' how the model matrix, X of the formula in \code{\link{DESeq}}, is
#' formed. "standard" is as created by \code{model.matrix} using the
#' design formula. "expanded" includes an indicator variable for each
#' level of factors with 3 or more levels, in order to ensure symmetric
#' behavior of the prior on log2 fold changes. betaPrior must be set
#' to TRUE in order for expanded model matrices to be fit.
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
nbinomWaldTest <- function(object, betaPrior=TRUE, betaPriorVar, modelMatrixType,
                           maxit=100, useOptim=TRUE, quiet=FALSE,
                           useT=FALSE, df, useQR=TRUE) {
  if (is.null(dispersions(object))) {
    stop("testing requires dispersion estimates, first call estimateDispersions()")
  }
  stopifnot(is.logical(betaPrior))
  if (missing(modelMatrixType)) {
    if (factorPresentThreeOrMoreLevels(object) & betaPrior) {
      modelMatrixType <- "expanded"
    } else {
      modelMatrixType <- "standard"
    }
  }
  if (modelMatrixType == "expanded" & !betaPrior) {
    stop("expanded model matrices require a beta prior")
  }
  stopifnot(length(maxit)==1)
  if ("results" %in% mcols(mcols(object))$type) {
    if (!quiet) message("you had results columns, replacing these")
    object <- removeResults(object)
  }
  if (!"allZero" %in% names(mcols(object))) {
    object <- getBaseMeansAndVariances(object)
  }
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]

  # if there are interaction terms present in the design
  # then we should only use the prior on the interaction terms
  termsOrder <- attr(terms.formula(design(object)),"order")
  interactionPresent <- any(termsOrder > 1)
  if (any(termsOrder > 2) & modelMatrixType == "expanded") {
    stop("interactions higher than 2nd order and usage of expanded model matrices
has not been implemented. we recommend instead using a likelihood
ratio test, i.e. DESeq with argument test='LRT'.")
  }
  priorOnlyInteraction <- interactionPresent & betaPrior & missing(betaPriorVar)

  if (!betaPrior) {
    # fit the negative binomial GLM without a prior
    # (in actuality a very wide prior with standard deviation 1e3 on log2 fold changes)
    fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim, useQR=useQR)
    H <- fit$hat_diagonals
    modelMatrix <- fit$modelMatrix
    modelMatrixNames <- fit$modelMatrixNames
    # record the wide prior variance which was used in fitting
    betaPriorVar <- rep(1e6, ncol(fit$modelMatrix))
  }

  if (betaPrior) {
    # first, fit the negative binomial GLM without a prior,
    # used to construct the prior variances
    # and for the hat matrix diagonals for calculating Cook's distance
    if (modelMatrixType == "standard") {
      fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim, useQR=useQR)
      modelMatrix <- fit$modelMatrix
      modelMatrixNames <- fit$modelMatrixNames
    } else {
      # expanded model matrices: want to make sure the prior
      # doesn't change from re-leveling
      modelMatrix <- makeReleveledModelMatrix(object)
      modelMatrixNames <- colnames(modelMatrix)
      fit <- fitNbinomGLMs(objectNZ, modelMatrix=modelMatrix,
                           maxit=maxit, useOptim=useOptim,
                           useQR=useQR, renameCols=FALSE)
    }
    H <- fit$hat_diagonal
    betaMatrix <- fit$betaMatrix
    colnames(betaMatrix) <- modelMatrixNames
    
    if (missing(betaPriorVar)) {
      # estimate the variance of the prior on betas
      # if expanded, first calculate LFC for all possible contrasts
      if (modelMatrixType == "expanded") {
        betaMatrix <- addAllContrasts(object, betaMatrix)
      } 
      if (nrow(fit$betaMatrix) > 1) {
        betaPriorVar <- apply(betaMatrix, 2, function(x) {
          # infinite betas are halted when |beta| > 10
          # so this test removes them
          useSmall <- abs(x) < 8
          # if no more betas pass test, return wide prior
          if (sum(useSmall) == 0 ) {
            return(1e6)
          } else {
            mean(x[useSmall]^2)
          }
        }) 
      } else {
        betaPriorVar <- (betaMatrix)^2
      }
      names(betaPriorVar) <- colnames(betaMatrix)

      # find the names of betaPriorVar which correspond
      # to non-interaction terms and set these to a wide prior
      if (priorOnlyInteraction) {
        nonInteractionCols <- getNonInteractionColumnIndices(object, modelMatrix)
        if (modelMatrixType == "standard") widePrior <- 1e6 else widePrior <- 1e3
        betaPriorVar[nonInteractionCols] <- widePrior
        if (modelMatrixType == "expanded") {
          # also set a wide prior for additional contrasts which were added
          # for calculation of the prior variance in the case of
          # expanded model matrices
          designFactors <- getDesignFactors(object)
          betaPriorVar[names(betaPriorVar) %in% paste0(designFactors,"Cntrst")] <- widePrior
        }
      }
      # intercept set to wide prior
      if ("Intercept" %in% modelMatrixNames) {
        betaPriorVar[which(names(betaPriorVar) == "Intercept")] <- 1e6
      }
      
      if (modelMatrixType == "expanded") {
        # bring over beta priors from the GLM fit without prior.
        # for factors: prior variance of each level are the average of the
        # prior variances for the levels present in the previous GLM fit
        betaPriorExpanded <- averagePriorsOverLevels(object, betaPriorVar)
        betaPriorVar <- betaPriorExpanded
      }
    } else {
      # else we are provided the prior variance:
      # check if the lambda is the correct length
      # given the design formula
      if (modelMatrixType == "expanded") {
        modelMatrix <- makeExpandedModelMatrix(object)
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
      modelMatrixNames <- fit$modelMatrixNames
    } else {
      modelMatrix <- makeExpandedModelMatrix(object)
      modelMatrixNames <- colnames(modelMatrix)
      fit <- fitNbinomGLMs(objectNZ, lambda=lambda, maxit=maxit, useOptim=useOptim,
                           useQR=useQR, modelMatrix=modelMatrix, renameCols=FALSE)
    }
  }

  # store mu in case the user did not call estimateDispersionsGeneEst
  assays(objectNZ)[["mu"]] <- fit$mu
  assays(object)[["mu"]] <- buildMatrixWithNARows(fit$mu, mcols(object)$allZero)

  # store the prior variance directly as an attribute
  # of the DESeqDataSet object, so it can be pulled later by
  # the results function (necessary for setting max Cook's distance)
  attr(object,"betaPrior") <- betaPrior
  attr(object,"betaPriorVar") <- betaPriorVar
  attr(object,"modelMatrix") <- modelMatrix
  attr(object,"modelMatrixType") <- modelMatrixType

  m <- nrow(modelMatrix)
  p <- ncol(modelMatrix)

  # calculate Cook's distance
  cooks <- calculateCooksDistance(objectNZ, H, p)

  # record maximum Cook's
  maxCooks <- recordMaxCooks(design(object), colData(object), fit$modelMatrix, cooks, nrow(objectNZ))

  # store Cook's distance for each sample
  assays(object)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)
  
  # add betas, standard errors and Wald p-values to the object
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
  
  resultsList <- c(matrixToList(betaMatrix),
                   matrixToList(betaSE),
                   matrixToList(WaldStatistic),
                   matrixToList(WaldPvalue),
                   list(betaConv = betaConv,
                        betaIter = fit$betaIter,
                        deviance = -2 * fit$logLike,
                        maxCooks = maxCooks))
  WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
  
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  if (betaPrior) {
    coefInfo <- paste("log2 fold change (MAP):",modelMatrixNamesSpaces)
  } else {
    coefInfo <- paste("log2 fold change:",modelMatrixNamesSpaces)
  }
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)
  statInfo <- paste("Wald statistic:",modelMatrixNamesSpaces)
  pvalInfo <- paste("Wald test p-value:",modelMatrixNamesSpaces)
  
  mcols(WaldResults) <- DataFrame(type = rep("results",ncol(WaldResults)),
                                  description = c(coefInfo, seInfo, statInfo, pvalInfo,
                                    "convergence of betas",
                                    "iterations for betas",
                                    "deviance for the fitted model",
                                    "maximum Cook's distance for row"))
  
  mcols(object) <- cbind(mcols(object),WaldResults)
  return(object)
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
#' \code{design(object)}
#' @param reduced a reduced formula to compare against, e.g.
#' the full model with a term or terms of interest removed
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
                      maxit=100, useOptim=TRUE, quiet=FALSE,
                      useQR=TRUE) {
  if (is.null(dispersions(object))) {
    stop("testing requires dispersion estimates, first call estimateDispersions()")
  }
  if (missing(reduced)) {
    stop("please provide a reduced formula for the likelihood ratio test, e.g. nbinomLRT(object, reduced = ~ 1)")
  }
  if (any(mcols(mcols(object))$type == "results")) {
    if (!quiet) message("you had results columns, replacing these")
    object <- removeResults(object)
  } 

  # try to form model matrices, test for difference
  # in residual degrees of freedom
  fullModelMatrix <- model.matrix(full,
                                  data=as.data.frame(colData(object)))
  reducedModelMatrix <- model.matrix(reduced,
                                     data=as.data.frame(colData(object)))
  df <- ncol(fullModelMatrix) - ncol(reducedModelMatrix)
  if (df < 1) stop("less than one degree of freedom, perhaps full and reduced models are not in the correct order")

  if (!"allZero" %in% names(mcols(object))) {
    object <- getBaseMeansAndVariances(object)
  }
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]
  
  fullModel <- fitNbinomGLMs(objectNZ, modelFormula=full, maxit=maxit,
                             useOptim=useOptim, useQR=useQR)
  reducedModel <- fitNbinomGLMs(objectNZ, modelFormula=reduced, maxit=maxit,
                                useOptim=useOptim, useQR=useQR)

  attr(object, "betaPrior") <- FALSE
  attr(object,"modelMatrix") <- fullModelMatrix
  attr(object,"modelMatrixType") <- "standard"
  
  p <- ncol(fullModelMatrix)
  m <- nrow(fullModelMatrix)
  H <- fullModel$hat_diagonals

  # store mu in case the user did not call estimateDispersionsGeneEst
  assays(objectNZ)[["mu"]] <- fullModel$mu
  assays(object)[["mu"]] <- buildMatrixWithNARows(fullModel$mu, mcols(object)$allZero)
  
  # calculate Cook's distance
  cooks <- calculateCooksDistance(objectNZ, H, p)

  # record maximum of Cook's
  maxCooks <- recordMaxCooks(design(object), colData(object), fullModelMatrix, cooks, nrow(objectNZ))

  # store Cook's distance for each sample
  assays(object)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)
  
  if (any(!fullModel$betaConv)) {
    if (!quiet) message(paste(sum(!fullModel$betaConv),"rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT"))
  }

  # calculate LRT statistic and p-values
  LRTStatistic <- (2 * (fullModel$logLike - reducedModel$logLike))
  LRTPvalue <- pchisq(LRTStatistic, df=df, lower.tail=FALSE)

  # continue storing LRT results
  resultsList <- c(matrixToList(fullModel$betaMatrix),
                   matrixToList(fullModel$betaSE),
                   list(LRTStatistic = LRTStatistic,
                        LRTPvalue = LRTPvalue,
                        fullBetaConv = fullModel$betaConv,
                        reducedBetaConv = reducedModel$betaConv,
                        betaIter = fullModel$betaIter,
                        deviance = -2 * fullModel$logLike,
                        maxCooks = maxCooks))
  LRTResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)

  modelComparison <- paste0("'",paste(as.character(full),collapse=" "),
                            "' vs '", paste(as.character(reduced),collapse=" "),"'")

  modelMatrixNames <- colnames(fullModel$betaMatrix)
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  coefInfo <- paste("log2 fold change:",modelMatrixNamesSpaces)
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)
  statInfo <- paste("LRT statistic:",modelComparison)
  pvalInfo <- paste("LRT p-value:",modelComparison)

  mcols(LRTResults) <- DataFrame(type = rep("results",ncol(LRTResults)),
                                 description = c(coefInfo, seInfo,
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
#' of replicates (including that sample). For more information on Cook's distance, please see the two
#' sections of the vignette: 'Dealing with count outliers' and 'Count outlier detection'.
#' 
#' @param dds a DESeqDataSet object, which has already been processed by
#' either DESeq, nbinomWaldTest or nbinomLRT, and therefore contains a matrix
#' 'cooks' contained in assays(dds). These are the Cook's distances which will
#' be used to define outlier counts.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from
#' each end of the normalized counts for a gene before the mean is computed
#' @param cooksCutoff the threshold for defining an outlier to be replaced.
#' Defaults to the .99 quantile of the F(p, m - p) distribution, where p is
#' the number of parameters and m is the number of samples.
#' @param minReplicates the number of replicate samples necessary to consider
#' a sample eligible for replacement. Outlier counts will not be replaced
#' if the sample is in a cell which has less than minReplicates replicates.
#' @param whichSamples a numeric or logical index of which samples will
#' have outliers replaced, if missing, minReplicates is used.
#'
#' @seealso \code{\link{DESeq}}
#'
#' @return a DESeqDataSet with replaced counts in the slot returned by
#' \code{\link{counts}} and the original counts preserved in
#' \code{assays(dds)[["originalCounts"]]}
#' 
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(n=100)
#' dds <- DESeq(dds)
#' ddsReplace <- replaceOutliersWithTrimmedMean(dds)
#'
#' @export
replaceOutliersWithTrimmedMean <- function(dds,trim=.2,cooksCutoff,minReplicates=7,whichSamples) {
  if (is.null(attr(dds,"modelMatrix")) | !("cooks" %in% names(assays(dds)))) {
    stop("first run DESeq, nbinomWaldTest, or nbinomLRT to identify outliers")
  }
  if (minReplicates < 3) {
    stop("at least 3 replicates are necessary in order to indentify a sample as a count outlier")
  }
  p <- ncol(attr(dds,"modelMatrix"))
  m <- ncol(dds)
  if (m <= p) {
    assays(dds)[["originalCounts"]] <- counts(dds)
    return(dds)
  }
  if (missing(cooksCutoff)) {
    cooksCutoff <- qf(.99, p, m - p)
  }
  idx <- which(assays(dds)[["cooks"]] > cooksCutoff)
  trimBaseMean <- apply(counts(dds,normalized=TRUE),1,mean,trim=trim)
  # build a matrix of counts based on the trimmed mean and the size factors
  if (!is.null(normalizationFactors(dds))) {
    replacementCounts <- as.integer(matrix(rep(trimBaseMean,ncol(dds)),ncol=dds) * 
                                    normalizationFactors(dds))
  } else {
    replacementCounts <- as.integer(outer(trimBaseMean,
                                          sizeFactors(dds), "*"))
  }
  # replace only those values which fall above the cutoff on Cook's distance
  newCounts <- counts(dds)
  newCounts[idx] <- replacementCounts[idx]
  if (missing(whichSamples)) {
    whichSamples <- nOrMoreInCell(attr(dds,"modelMatrix"), n = minReplicates)
  }
  if (is.logical(whichSamples)) whichSamples <- which(whichSamples)
  assays(dds)[["originalCounts"]] <- counts(dds)
  if (length(whichSamples) == 0) {
    return(dds)
  }
  counts(dds)[,whichSamples] <- newCounts[,whichSamples]
  dds
}



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
                                  description = c("the base mean over all rows",
                                    "the base variance over all rows",
                                    "all counts in a row are zero"))
  if (all(c("baseMean","baseVar","allZero") %in% names(mcols(object)))) {
      mcols(object)[c("baseMean","baseVar","allZero")] <- meanVarZero
  } else {
      mcols(object) <- cbind(mcols(object),meanVarZero)
  }
  return(object)
}


# Estimate a parametric fit of dispersion to the mean intensity
parametricDispersionFit <- function( means, disps ) {
   coefs <- c( .1, 1 )
   iter <- 0
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which( (residuals > 1e-4) & (residuals < 15) )
      fit <- glm( disps[good] ~ I(1/means[good]),
         family=Gamma(link="identity"), start=coefs )
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      if( !all( coefs > 0 ) )
         stop(simpleError("parametric dispersion fit failed"))
      if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
         break
      iter <- iter + 1
      if( iter > 10 ) {
         stop(simpleError("dispersion fit did not converge"))
         break }
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
  fit <- locfit(logDisps ~ logMeans, data=d[disps >= minDisp*10,],
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
#
# return a list of results, with coefficients and standard
# errors on the log2 scale
fitNbinomGLMs <- function(object, modelMatrix, modelFormula, alpha_hat, lambda,
                          renameCols=TRUE, betaTol=1e-8, maxit=100, useOptim=TRUE,
                          useQR=TRUE, forceOptim=FALSE) {
  if (missing(modelFormula)) {
    modelFormula <- design(object)
  }
  if (missing(modelMatrix)) {
   modelMatrix <- model.matrix(modelFormula, data=as.data.frame(colData(object)))
  }
  modelMatrixNames <- colnames(modelMatrix)

  if (renameCols) {
    convertNames <- renameModelMatrixColumns(modelMatrixNames,
                                             as.data.frame(colData(object)),
                                             modelFormula)
    convertNames <- convertNames[convertNames$from %in% modelMatrixNames,]
    modelMatrixNames[match(convertNames$from, modelMatrixNames)] <- convertNames$to
  }
  
  modelMatrixNames[modelMatrixNames == "(Intercept)"] <- "Intercept"
  
  if (!is.null(normalizationFactors(object))) {
    normalizationFactors <- normalizationFactors(object)
  } else { 
    normalizationFactors <- matrix(rep(sizeFactors(object),each=nrow(object)),
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
  if (modelFormula == formula(~ 1) & all(lambda <= 1e-6)) {
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
                  modelMatrix=modelMatrix, modelMatrixNames = modelMatrixNames,
                  nterms=1, hat_diagonals=hat_diagonals)
      return(res)
  }
  
  if ("Intercept" %in% modelMatrixNames) {
    beta_mat <- matrix(0, ncol=ncol(modelMatrix), nrow=nrow(object))
    beta_mat[,which(modelMatrixNames == "Intercept")] <- log(mcols(object)$baseMean)
  } else {
    beta_mat <- matrix(1, ncol=ncol(modelMatrix), nrow=nrow(object))
  }
  
  # here we convert from the log2 scale of the betas
  # and the beta prior variance to the log scale
  # used in fitBeta.
  # so we divide by the square of the
  # conversion factor, log(2)
  lambdaLogScale <- lambda / log(2)^2

  betaRes <- fitBeta(ySEXP = counts(object), xSEXP = modelMatrix,
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
  betaSE <- log2(exp(1))*sqrt(betaRes$beta_var_mat)
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)

  # switch based on whether we should also use optim
  # on rows which did not converge
  if (useOptim) {
    rowsForOptim <- which(!betaConv | !rowStable | !rowVarPositive)
  } else {
    rowsForOptim <- which(!rowStable | !rowVarPositive)
  }
  if (forceOptim) {
    rowsForOptim <- seq_along(betaConv)
  }
  
  if (length(rowsForOptim) > 0) {
    scaleCols <- apply(modelMatrix,2,function(z) max(abs(z)))
    x <- sweep(modelMatrix,2,scaleCols,"/")
    lambdaColScale <- lambda / scaleCols^2
    lambdaColScale <- ifelse(lambdaColScale == 0, 1e-6, lambdaColScale)
    lambdaLogScaleColScale <- lambdaLogScale / scaleCols^2
    large <- 30
    for (row in rowsForOptim) {
      if (rowStable[row]) {
        betaRow <- betaMatrix[row,] * scaleCols
      } else {
        betaRow <- beta_mat[row,] * scaleCols
      }
      betaRow <- pmin(pmax(betaRow, -large), large)
      nf <- normalizationFactors[row,]
      k <- counts(object)[row,]
      alpha <- alpha_hat[row]
      objectiveFn <- function(p) {
        mu_row <- as.numeric(nf * 2^(x %*% p))
        prior <- sum(dnorm(p,0,sqrt(1/lambdaColScale),log=TRUE))
        logLike <- sum(dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE))
        softBox <- sum( ifelse(abs(p) < 30, 0, (abs(p) - 30)) )
        # -1 times the posterior plus the soft box penalty
        -1 * (logLike + prior) + softBox
      }
      o <- optim(betaRow, objectiveFn, method="L-BFGS-B")
      if (length(lambdaLogScale) > 1) {
        ridge <- diag(lambdaLogScaleColScale)
      } else {
        ridge <- as.matrix(lambdaLogScaleColScale,ncol=1)
      }
      # if we converged, change betaConv to TRUE
      if (o$convergence == 0) {
        betaConv[row] <- TRUE
      }
      # with or without convergence, store the estimate from optim
      betaMatrix[row,] <- o$par / scaleCols
      # calculate the standard errors
      mu_row <- as.numeric(nf * 2^(x %*% o$par))
      w <- diag((mu_row^-1 + alpha)^-1)
      xtwx <- t(x) %*% w %*% x
      xtwxRidgeInv <- solve(xtwx + ridge)
      sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
      betaSE[row,] <- log2(exp(1)) * sqrt(diag(sigma)) / scaleCols
      # store the new mu vector
      mu[row,] <- mu_row
      logLike[row] <- sum(dnbinom(k, mu=mu_row, size=1/alpha, log=TRUE))
    }
  }
  
  list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter,
       deviance = betaRes$deviance,
       modelMatrix=modelMatrix, modelMatrixNames = modelMatrixNames,
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
# kappa_0SEXP a parameter used in calculting the initial proposal
#   for the backtracking search
#   initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
# use_priorSEXP boolean variable, whether to use a prior or just calculate the MLE
#
# return a list with elements: log_alpha, iter, iter_accept, last_change, initial_lp, intial_dlp, last_lp, last_dlp, last_d2lp
fitDisp <- function (ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP, log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP, tolSEXP, maxitSEXP, use_priorSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitDisp))
  na.test <- sapply(list(ySEXP, xSEXP, mu_hatSEXP, log_alphaSEXP, log_alpha_prior_meanSEXP, log_alpha_prior_sigmasqSEXP, min_log_alphaSEXP, kappa_0SEXP, tolSEXP, maxitSEXP, use_priorSEXP), function(x) any(is.na(x)))
  if (any(na.test)) {
    stop(paste("in call to fitDisp, the following arguments contain NA:",paste(arg.names[na.test],collapse=", ")))
  }
  .Call("fitDisp", ySEXP=ySEXP, xSEXP=xSEXP, mu_hatSEXP=mu_hatSEXP,
        log_alphaSEXP=log_alphaSEXP, log_alpha_prior_meanSEXP=log_alpha_prior_meanSEXP,
        log_alpha_prior_sigmasqSEXP=log_alpha_prior_sigmasqSEXP,
        min_log_alphaSEXP=min_log_alphaSEXP, kappa_0SEXP=kappa_0SEXP,
        tolSEXP=tolSEXP, maxitSEXP=maxitSEXP, use_priorSEXP=use_priorSEXP,
        PACKAGE = "DESeq2")
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
fitBeta <- function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP, beta_matSEXP, lambdaSEXP, tolSEXP, maxitSEXP, useQRSEXP) {
  if ( missing(contrastSEXP) ) {
    contrastSEXP <- c(1,rep(0,ncol(xSEXP)-1))
  }
  # test for any NAs in arguments
  arg.names <- names(formals(fitBeta))
  na.test <- sapply(list(ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, contrastSEXP, beta_matSEXP, lambdaSEXP, tolSEXP, maxitSEXP, useQRSEXP), function(x) any(is.na(x)))
  if (any(na.test)) {
    stop(paste("in call to fitBeta, the following arguments contain NA:",paste(arg.names[na.test],collapse=", ")))
  }
  .Call("fitBeta", ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP,
        alpha_hatSEXP=alpha_hatSEXP, contrastSEXP=contrastSEXP,
        beta_matSEXP=beta_matSEXP,
        lambdaSEXP=lambdaSEXP, tolSEXP=tolSEXP, maxitSEXP=maxitSEXP,
        useQRSEXP=useQRSEXP,
        PACKAGE = "DESeq2")
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


# calculate a robust method of moments dispersion
# using the squared MAD and the row mean
# the idea is to get a sense of the dispersion
# excluding individual outlier counts which would raise the estimate
robustMethodOfMomentsDisp <- function(counts) {
  v <- apply(counts,1,mad)^2
  m <- rowMeans(counts)
  alpha <- pmax(v - m, 0) / m^2
  alpha <- ifelse(alpha > 0, alpha, min(alpha))
}


calculateCooksDistance <- function(object, H, p) {
  dispersions <- robustMethodOfMomentsDisp(counts(object,normalized=TRUE))
  V <- assays(object)[["mu"]] + dispersions * assays(object)[["mu"]]^2
  PearsonResSq <- (counts(object) - assays(object)[["mu"]])^2 / V
  cooks <- PearsonResSq / p  * H / (1 - H)^2
  cooks
}


# this function breaks out the logic for calculating the max Cook's distance:
# the samples over which max Cook's distance is calculated:
#
# if all the variables in the design are factors, then those samples with 3 or more replicates per cell
# if one or more are not factor, then those samples such that the matrix is full rank after removing the row
#
# if m == p or there are no samples over which to calculate max Cook's, then give NA
recordMaxCooks <- function(design, colData, modelMatrix, cooks, numRow) {
    if (allFactors(design, colData)) {
        samplesForCooks <- nOrMoreInCell(modelMatrix, n=3)
    } else {
        samplesForCooks <- leaveOneOutFullRank(modelMatrix)
    }
    p <- ncol(modelMatrix)
    m <- nrow(modelMatrix)
    if ((m > p) & any(samplesForCooks)) {
        maxCooks <- apply(cooks[,samplesForCooks,drop=FALSE], 1, max)
    } else {
        maxCooks <- rep(NA, numRow)
    }
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

# are all the variables in the design matrix factors?
allFactors <- function(design, colData) {
    designVars <- all.vars(formula(design))
    designVarsClass <- sapply(designVars, function(v) class(colData[[v]]))
    all(designVarsClass == "factor")
}

# returns TRUE or FALSE if removing row would leave matrix full rank
leaveOneOutFullRank <- function(modelMatrix) {
  sapply(seq_len(nrow(modelMatrix)), function(i) qr(modelMatrix[-i,])$rank) == ncol(modelMatrix)
}

propMaxOfTotal <- function(counts, pseudocount=8) {
    (apply(counts, 1, max) + pseudocount)/(rowSums(counts) + pseudocount*ncol(counts))
}



# an unexported diagnostic function
# to retrieve the covariance matrix
# for the GLM coefficients of a single row
# only for standard model matrices
covarianceMatrix <- function(object, rowNumber) {
  # convert coefficients to log scale
  coefColumns <- names(mcols(object))[grep("log2 fold change",mcols(mcols(object))$description)]
  beta <- log(2) * as.numeric(as.data.frame(mcols(object)[rowNumber,coefColumns]))
  x <- model.matrix(design(object), as.data.frame(colData(object)))
  y <- counts(object)[rowNumber,]
  sf <- sizeFactors(object)
  alpha <- dispersions(object)[rowNumber]
  mu.hat <- as.vector(sf * exp(x %*% beta))
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
  designVars <- all.vars(formula(design))
  designVarsClass <- sapply(designVars, function(v) class(colData(object)[[v]]))
  designVars[designVarsClass == "factor"]
}

factorPresentThreeOrMoreLevels <- function(object) {
  designFactors <- getDesignFactors(object)
  threeOrMore <- sapply(designFactors,function(v) nlevels(colData(object)[[v]]) >= 3)
  any(threeOrMore)
}
