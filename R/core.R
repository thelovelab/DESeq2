#' Differential expression analysis based on the negative binomial distribution
#'
#' This function performs a default analysis by calling, in order, the functions:
#' \code{\link{estimateSizeFactors}},
#' \code{\link{estimateDispersions}},
#' \code{\link{nbinomWaldTest}}.
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
#' column of the model matrix \eqn{X}{X}.  The sample-specific size factors
#' can be replaced by gene-specific normalization factors for each sample using
#' \code{\link{normalizationFactors}}.  For details on the fitting of the log2
#' fold changes and calculation of p-values see \code{\link{nbinomWaldTest}}.
#'
#' @return a \code{\link{DESeqDataSet}} object with results stored as
#' metadata columns.  The results can be accessed by calling the \code{\link{results}}
#' function.  By default this will return the log2 fold changes and p-values for the last
#' variable in the design formula.  See \code{\link{results}} for how to access results
#' for other variables.
#'
#' @param object a DESeqDataSet object, see the constructor functions
#' \code{\link{DESeqDataSet}},
#' \code{\link{DESeqDataSetFromMatrix}},
#' \code{\link{DESeqDataSetFromHTSeqCount}}.
#' @param fitType either "parametric", "local", or "mean"
#' for the type of fitting of dispersions to the mean intensity.
#' See \code{\link{estimateDispersions}} for description.
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients (Tikhonov/ridge regularization)
#' See \code{\link{nbinomWaldTest}} for description.
#' @param pAdjustMethod the method to use for adjusting p-values, see \code{?p.adjust}
#' @param quiet whether to print messages at each step
#'
#' @author Michael Love
#'
#' @references Simon Anders, Wolfgang Huber: Differential expression analysis for sequence count data. Genome Biology 11 (2010) R106, \url{http://dx.doi.org/10.1186/gb-2010-11-10-r106}
#'
#' @import Biobase GenomicRanges IRanges
#' @importFrom locfit locfit
#' @importFrom genefilter rowVars
#' @importFrom RColorBrewer brewer.pal
#' @useDynLib DESeq2
#'
#' @seealso \code{\link{nbinomWaldTest}}, \code{\link{nbinomLRT}}
#'
#' @examples
#'
#' dds <- makeExampleDESeqDataSet(betaSd=1)
#' dds <- DESeq(dds)
#' res <- results(dds)
#'
#' @export
DESeq <- function(object,fitType=c("parametric","local","mean"),betaPrior=TRUE,pAdjustMethod="BH",quiet=FALSE) {
  if (!quiet) message("estimating size factors")
  object <- estimateSizeFactors(object)
  if (!quiet) message("estimating dispersions")
  object <- estimateDispersions(object,fitType=fitType, quiet=quiet)
  if (!quiet) message("fitting generalized linear model")
  object <- nbinomWaldTest(object, betaPrior=betaPrior,
                           pAdjustMethod=pAdjustMethod, quiet=quiet)
  object
}

#' Make a simulated DESeqDataSet
#'
#' Constructs a simulated dataset of negative binomial data from
#' two conditions. By default, there are no fold changes between
#' the two conditions, but this can be adjusted with the \code{betaSd} argument.
#'
#' @param n number of rows
#' @param m number of columns
#' @param betaSd the standard deviation for non-intercept betas, i.e. beta ~ N(0,betaSd)
#' @param interceptMean the mean of the intercept betas (log2 scale)
#' @param interceptSd the standard deviation of the intercept betas (log2 scale)
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
makeExampleDESeqDataSet <- function(n=1000,m=10,betaSd=0,interceptMean=4,interceptSd=2,dispMeanRel=function(x) 4/x + .1,sizeFactors=rep(1,m)) {
  beta <- cbind(rnorm(n,interceptMean,interceptSd),rnorm(n,0,betaSd))
  dispersion <- dispMeanRel(2^(beta[,1]))
  colData <- DataFrame(sample = paste("sample",1:m,sep=""),
                       condition=factor(rep(c("A","B"),times=c(ceiling(m/2),floor(m/2)))))
  x <- model.matrix(~ colData$condition)
  mu <- 2^(t(x %*% t(beta))) * rep(sizeFactors, each=n)
  countData <- matrix(rnbinom(m*n, mu=mu, size=1/dispersion), ncol=m)
  rownames(colData) <- colData$sample
  rowData <- GRanges("1",IRanges(start=(1:n - 1) * 100 + 1,width=100))
  names(rowData) <- paste0("feature",1:n)
  object <- DESeqDataSetFromMatrix(countData = countData,
                                                colData = colData,
                                                design = formula(~ condition),
                                                rowData = rowData)
  trueVals <- DataFrame(trueIntercept = beta[,1],
                        trueBeta = beta[,2],
                        trueDisp = dispersion)
  mcols(trueVals) <- DataFrame(type=rep("input",ncol(trueVals)),
                               description=c("simulated intercept values",
                                 "simulated beta values",
                                 "simulated dispersion values"))
  mcols(object) <- cbind(mcols(object),trueVals)
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
#' @return a vector with the estimates size factors, one element per column
#' @author Simon Anders
#' @seealso \code{\link{estimateSizeFactors}}
#' @examples
#' 
#' dds <- makeExampleDESeqDataSet()
#' estimateSizeFactorsForMatrix(counts(dds))
#' 
#' @export
estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log( counts ) )
   if (all(is.infinite(loggeomeans))) {
     stop("every gene contains at least one zero, cannot compute log geometric means")
   }
   apply( counts, 2, function(cnts)
      exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
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
#' @param priorVar the variance of the normal prior on the log dispersions.
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
#' mcols(dds)$dispFit <- rep(median(mcols(dds)$dispGeneEst,na.rm=TRUE),nrow(dds))
#' dds <- estimateDispersionsMAP(dds)
#' plotDispEsts(dds)
#'
#' @seealso \code{\link{estimateDispersions}}
#'
#' @export
estimateDispersionsGeneEst <- function(object, minDisp=1e-8, kappa_0=1, dispTol=1e-6, maxit=100, quiet=FALSE) {
  if ("dispGeneEst" %in% names(mcols(object))) {
    if (!quiet) message("you had estimated gene-wise dispersions, removing these")
    mcols(object) <- mcols(object)[,!names(mcols(object))  == "dispGeneEst"]
  }
  if (log(minDisp/10) <= -30) {
    stop("for computational stability, log(minDisp/10) should be above -30")
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

  # fitNbinomGLMs returns mu and modelMatrix
  fit <- fitNbinomGLMs(objectNZ, alpha_hat=alpha_hat)
  
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
estimateDispersionsFit <- function(object,fitType=c("parametric","local","mean"),minDisp=1e-8, quiet=FALSE) {
  if ("dispFit" %in% names(mcols(object))) {
    if (!quiet) message("you had estimated fitted dispersions, removing these")
    mcols(object) <- mcols(object)[,!names(mcols(object)) == "dispFit"]
  }
  objectNZ <- object[!mcols(object)$allZero,]
  useForFit <- mcols(objectNZ)$dispGeneEstConv

  # take the first fitType
  fitType <- fitType[1]
  
  if (fitType == "parametric") {
    trial <- try(dispFunction <- parametricDispersionFit(mcols(objectNZ)$baseMean[useForFit],
                                                         mcols(objectNZ)$dispGeneEst[useForFit]),
                 silent=TRUE)
    if (!inherits(trial,"try-error")) {
      dispFit <- dispFunction(mcols(objectNZ)$baseMean)
    } else {
      warning("parametric fit failed, trying local fit. use plotDispEsts() to check the quality of the local fit")
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
estimateDispersionsMAP <- function(object, outlierSD=2, priorVar, minDisp=1e-8, kappa_0=1, dispTol=1e-6, maxit=100, quiet=FALSE) {
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
    priorVarCalc <- pmax((varLogDispEsts - expVarLogDisp), .25)
  } else {
    # we have m = p, so do not try to substract sampling variance
    priorVarCalc <- varLogDispEsts
  }
  

  # fill in the calculated prior variance
  if (missing(priorVar)) {
    priorVar <- priorVarCalc
  }
  attr( dispersionFunction(object), "priorVar" ) <- priorVar
  
  # set prior variance for fitting dispersion
  log_alpha_prior_sigmasq <- priorVar

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

  if (any(!resultsList$dispConv)) {
    if (!quiet) message(paste(sum(!resultsList$dispConv),"rows did not converge in dispersion, labelled in mcols(object)$dispConv. Use larger maxit argument with estimateDispersions"))
  }
  
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
#' non-intercept coefficient is calculated as the mean squared maximum likelihood
#' estimates over the genes which do not contain zeros for some condition;
#' the final coefficients are then maximum a posteriori estimates
#' (using Tikhonov/ridge regularization) using this prior.
#' The use of a prior has little effect on genes with high counts and helps to
#' moderate the large spread in coefficients for genes with low counts.
#'
#' For calculating Wald test p-values, the coefficients are scaled by their
#' standard errors and then compared to a normal distribution.  From
#' examination of Wald statistics for real datasets, the effect of the
#' prior on dispersion estimates results in a Wald statistic
#' distribution which is approximately normal.
#'
#' Cook's distance for each sample are accessible as a matrix "cooks"
#' stored in the assays() list. This measure is useful for identifying rows where the
#' observed counts might not fit to a negative binomial distribution.
#' 
#' The Wald test can be replaced with the \code{\link{nbinomLRT}}
#' for an alternative test of significance.
#'
#' @param object a DESeqDataSet
#' @param betaPrior whether or not to put a zero-mean normal prior on
#' the non-intercept coefficients (Tikhonov/ridge regularization)
#' @param pAdjustMethod the method to use for adjusting p-values, see \code{?p.adjust}
#' @param priorSigmaSq a vector with length equal to the number of
#' model terms including the intercept.
#  priorSigmaSq gives the variance of the prior on the sample betas,
#' which if missing is estimated from the rows which do not have any
#' zeros
#' @param cooksCutoff theshold on Cook's distance, such that if one or more
#' samples for a row have a distance higher, the p-value for the row is
#' set to NA.
#' The default cutoff is the .75 quantile of the F(p, m-p) distribution,
#' where p is the number of coefficients being fitted and m is the number of samples.
#' Set to Inf or FALSE to disable the resetting of p-values to NA.
#' Note: this test excludes the Cook's distance of samples whose removal
#' would result in rank deficient design matrix.
#' @param maxit the maximum number of iterations to allow for convergence of the
#' coefficient vector
#' @param useOptim whether to use the native optim function on rows which do not
#' converge within maxit
#' @param quiet whether to print messages at each step
#' @param useT whether to use a t-distribution as a null distribution,
#' for significance testing of the Wald statistics.
#' If FALSE, a standard normal null distribution is used.
#' @param df the degrees of freedom for the t-distribution
#'
#' @return a DESeqDataSet with results columns accessible
#' with the \code{\link{results}} function.  The coefficients and standard errors are
#' reported on a log2 scale.
#'
#' @seealso \code{\link{nbinomLRT}}
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
nbinomWaldTest <- function(object, betaPrior=TRUE, pAdjustMethod="BH", priorSigmaSq, cooksCutoff,
                           maxit=100, useOptim=TRUE, quiet=FALSE, useT=FALSE, df) {
  if (is.null(dispersions(object))) {
    stop("testing requires dispersion estimates, first call estimateDispersions()")
  }  
  if ("results" %in% mcols(mcols(object))$type) {
    if (!quiet) message("you had results columns, replacing these")
    object <- removeResults(object)
  }
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]
  
  if (!betaPrior) {
    fit <- fitNbinomGLMs(objectNZ, maxit=maxit, useOptim=useOptim)
    H <- fit$hat_diagonals
  }
  
  if (betaPrior) {
    # we need the MLE betas to fit the prior variance
    # and for the hat matrix diagonals in order to
    # calculate Cook's distance
    fit <- fitNbinomGLMs(objectNZ, maxit=maxit)
    H <- fit$hat_diagonals
    if (missing(priorSigmaSq)) {
      # estimate the width of the prior on betas
      # excluding those rows with coefficients going to infinity.
      # if all betas are large, use a very large prior
      priorSigmaSq <- apply(fit$betaMatrix, 2, function(x) {
        useSmall <- abs(x) < 8
        if (sum(useSmall) == 0 ) {
          return(1e6)
        } else {
          mean(x[useSmall]^2)
        }
      })
      # except for intercept which we set to wide prior
      if ("Intercept" %in% fit$modelMatrixNames) {
        priorSigmaSq[which(fit$modelMatrixNames == "Intercept")] <- 1e6
      }
    } else {
      # we are provided the prior variance:
      # check if the lambda is the correct length
      # given the design formula
      p <- ncol(fit$modelMatrix)
      if (length(priorSigmaSq) != p) {
        stop(paste0("priorSigmaSq should have length",p))
      }
    }
    lambda <- 1/priorSigmaSq
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda, maxit=maxit, useOptim=useOptim)
  }

  m <- nrow(fit$modelMatrix)
  p <- ncol(fit$modelMatrix)

  # calculate Cook's distance
  if (missing(cooksCutoff)) {
    cooksCutoff <- qf(.75, p, m - p)
  }
  if (is.logical(cooksCutoff) & cooksCutoff) {
    cooksCutoff <- qf(.75, p, m - p)
  }
  cooks <- calculateCooksDistance(objectNZ, H, p)
  looFullRank <- leaveOneOutFullRank(fit$modelMatrix)
  if (m > p) {
    maxCooks <- apply(cooks[,looFullRank,drop=FALSE], 1, max)
    cooksOutlier <- maxCooks > cooksCutoff
  } else {
    maxCooks <- rep(NA, nrow(objectNZ))
    cooksOutlier <- rep(FALSE, nrow(objectNZ))
  }

  # store Cook's distance for each sample
  assays(object)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)

  modelMatrixNames <- fit$modelMatrixNames
  
  # add betas, standard errors and Wald p-values to the object
  betaMatrix <- fit$betaMatrix
  colnames(betaMatrix) <- modelMatrixNames
  betaSE <- fit$betaSE
  colnames(betaSE) <- paste0("SE_",modelMatrixNames)
  WaldStatistic <- betaMatrix/betaSE
  colnames(WaldStatistic) <- paste0("WaldStatistic_",modelMatrixNames)
  
  # if useT is set to TRUE, use a t-distribution
  if (useT) {
    priorVar <- attr( dispersionFunction(object), "priorVar" )
    WaldPvalue <- 2*pt(abs(WaldStatistic),df=df,lower.tail=FALSE)
  } else {
    WaldPvalue <- 2*pnorm(abs(WaldStatistic),lower.tail=FALSE)
  }
  colnames(WaldPvalue) <- paste0("WaldPvalue_",modelMatrixNames)

  # Set to NA the p-values for genes that have one or more
  # samples with Cook's distance beyond the cutof
  if (is.numeric(cooksCutoff) | cooksCutoff) {
    WaldPvalue[cooksOutlier,] <- NA
  }
  
  # if more than 1 row, we adjust p-values
  if (nrow(WaldPvalue) > 1) {  
    WaldAdjPvalue <- apply(WaldPvalue,2,p.adjust,method=pAdjustMethod)
  } else {
    WaldAdjPvalue <- WaldPvalue
  }
  colnames(WaldAdjPvalue) <- paste0("WaldAdjPvalue_",modelMatrixNames)
  betaConv <- fit$betaConv

  if (any(!betaConv)) {
    if (!quiet) message(paste(sum(!betaConv),"rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest"))
  }
  
  resultsList <- c(matrixToList(betaMatrix),
                   matrixToList(betaSE),
                   matrixToList(WaldStatistic),
                   matrixToList(WaldPvalue),
                   matrixToList(WaldAdjPvalue),
                   list(betaConv = betaConv,
                        betaIter = fit$betaIter,
                        deviance = -2 * fit$logLike,
                        maxCooks = maxCooks,
                        cooksOutlier = cooksOutlier))
  WaldResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)
  
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  if (betaPrior) {
    coefInfo <- paste("log2 fold change (MAP):",modelMatrixNamesSpaces)
  } else {
    coefInfo <- paste("log2 fold change:",modelMatrixNamesSpaces)
  }
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)
  statInfo <- paste("Wald test:",modelMatrixNamesSpaces)
  pvalInfo <- paste("Wald test:",modelMatrixNamesSpaces)
  adjInfo <- paste(paste("Wald test,",pAdjustMethod,"adj.:"),
                   modelMatrixNamesSpaces)
  
  mcols(WaldResults) <- DataFrame(type = rep("results",ncol(WaldResults)),
                                  description = c(coefInfo, seInfo, statInfo, pvalInfo, adjInfo,
                                    "convergence of betas",
                                    "iterations for betas",
                                    "deviance for the fitted model",
                                    "maximum Cook's distance for row",
                                    "whether maxCooks > cooksCutoff"))
  
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
#' Cook's distance for each sample are accessible as a matrix "cooks"
#' stored in the assays() list. This measure is useful for identifying rows where the
#' observed counts might not fit to a negative binomial distribution.
#' 
#' @param object a DESeqDataSet
#' @param full the full model formula, this should be the formula in
#' \code{design(object)}
#' @param reduced a reduced formula to compare against, e.g.
#' the full model with a variable of interest removed
#' @param pAdjustMethod the method to use for adjusting p-values, see \code{?p.adjust}
#' @param cooksCutoff theshold on Cook's distance, such that if one or more
#' samples for a row have a distance higher, the p-value for the row is
#' set to NA.
#' The default cutoff is the .75 quantile of the F(p, m-p) distribution,
#' where p is the number of coefficients being fitted and m is the number of samples.
#' Set to Inf or FALSE to disable the resetting of p-values to NA.
#' Note: this test excludes the Cook's distance of samples whose removal
#' would result in rank deficient design matrix.
#' @param maxit the maximum number of iterations to allow for convergence of the
#' coefficient vector
#' @param useOptim whether to use the native optim function on rows which do not
#' converge within maxit
#' @param quiet whether to print messages at each step
#'
#' @return a DESeqDataSet with new results columns accessible
#' with the \code{\link{results}} function.  The coefficients and standard errors are
#' reported on a log2 scale.
#' 
#' @seealso \code{\link{nbinomWaldTest}}
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
nbinomLRT <- function(object, full=design(object), reduced, pAdjustMethod="BH", cooksCutoff,
                      maxit=100, useOptim=TRUE, quiet=FALSE ) {
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

  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]
  
  fullModel <- fitNbinomGLMs(objectNZ, modelFormula=full, maxit=maxit, useOptim=useOptim)
  reducedModel <- fitNbinomGLMs(objectNZ, modelFormula=reduced, maxit=maxit, useOptim=useOptim)

  p <- ncol(fullModelMatrix)
  m <- nrow(fullModelMatrix)
  H <- fullModel$hat_diagonals

  # calculate Cook's distance
  if (missing(cooksCutoff)) {
    cooksCutoff <- qf(.75, p, m - p)
  }
  if (is.logical(cooksCutoff) & cooksCutoff) {
    cooksCutoff <- qf(.75, p, m - p)
  }
  cooks <- calculateCooksDistance(objectNZ, H, p)
  looFullRank <- leaveOneOutFullRank(fullModelMatrix)
  if (m > p) {
    maxCooks <- apply(cooks[,looFullRank,drop=FALSE], 1, max)
    cooksOutlier <- maxCooks > cooksCutoff
  } else {
    maxCooks <- rep(NA, nrow(objectNZ))
    cooksOutlier <- rep(FALSE, nrow(objectNZ))
  }

  # store Cook's distance for each sample
  assays(object)[["cooks"]] <- buildMatrixWithNARows(cooks, mcols(object)$allZero)
  
  if (any(!fullModel$betaConv)) {
    if (!quiet) message(paste(sum(!fullModel$betaConv),"rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT"))
  }

  # calculate LRT statistic and p-values
  LRTStatistic <- (2 * (fullModel$logLike - reducedModel$logLike))
  LRTPvalue <- pchisq(LRTStatistic, df=df, lower.tail=FALSE)
  
  # Set to NA the p-values for genes that have one or more
  # samples with Cook's distance beyond the cutoff
  if (is.numeric(cooksCutoff) | cooksCutoff) {
    LRTPvalue[cooksOutlier] <- NA
  }

  # continue storing LRT results
  LRTAdjPvalue <- p.adjust(LRTPvalue,method=pAdjustMethod)
  resultsList <- c(matrixToList(fullModel$betaMatrix),
                   matrixToList(fullModel$betaSE),
                   list(LRTStatistic = LRTStatistic,
                        LRTPvalue = LRTPvalue,
                        LRTAdjPvalue = LRTAdjPvalue,
                        fullBetaConv = fullModel$betaConv,
                        reducedBetaConv = reducedModel$betaConv,
                        betaIter = fullModel$betaIter,
                        deviance = -2 * fullModel$logLike,
                        maxCooks = maxCooks,
                        cooksOutlier = cooksOutlier))
  LRTResults <- buildDataFrameWithNARows(resultsList, mcols(object)$allZero)

  modelComparison <- paste0("'",paste(as.character(full),collapse=" "), "' vs '", paste(as.character(reduced),collapse=" "),"'")

  modelMatrixNames <- colnames(fullModel$betaMatrix)
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  coefInfo <- paste("log2 fold change:",modelMatrixNamesSpaces)
  seInfo <- paste("standard error:",modelMatrixNamesSpaces)
  statInfo <- paste("LRT statistic:",modelComparison)
  pvalInfo <- paste("LRT p-value:",modelComparison)
  adjInfo <- paste("LRT p-value,",pAdjustMethod,"adj.")

  mcols(LRTResults) <- DataFrame(type = rep("results",ncol(LRTResults)),
                                 description = c(coefInfo, seInfo,
                                   statInfo, pvalInfo, adjInfo,
                                   "convergence of betas for full model",
                                   "convergence of betas for reduced model",
                                   "iterations for betas for full model",
                                   "deviance of the full model",
                                   "maximum Cook's distance for row",
                                   "whether maxCooks > cooksCutoff"))
  mcols(object) <- cbind(mcols(object),LRTResults)
  return(object)
}

#' Extract results from a DESeq analysis
#'
#' Extract results from a DESeq analysis giving base means across samples,
#' log2 fold changes, standard errors, p-values and adjusted p-values.
#' Find available names for results; Return an object with results columns
#' removed.
#'
#' Multiple results can be returned for an analysis with multifactor design,
#' so \code{results} takes an argument \code{name}
#'
#' The available names can be checked using \code{resultsNames};
#' these are combined variable names and factor levels, potentially 
#' with minor changes made by the \code{DataFrame} function on column names
#' (e.g. dashes into periods).
#'
#' By default, results for the last variable will be returned. Information on the variable
#' represented and the test used for p-values (Wald test or likelihood ratio test)
#' is stored in the metadata columns, accessible by calling \code{mcol}
#' on the object returned by \code{results}.
#'
#' Results can be removed from an object by calling \code{removeResults}
#'
#' @param object a DESeqDataSet, on which one
#' of the following functions has already been called:
#' \code{\link{DESeq}}, \code{\link{nbinomWaldTest}}, or \code{\link{nbinomLRT}}
#' @param name the name of the coefficient for which to report log2 fold changes,
#' p-values and FDR
#'
#' @return For \code{results}: a DataFrame of results columns with metadata
#' columns of coefficient and test information
#'
#' For \code{resultsNames}: the names of the columns available as results,
#' usually a combination of the variable name and a level
#'
#' For \code{removeResults}: the original object with results metadata columns removed
#'
#' @seealso \code{\link{DESeq}}
#'
#' @examples
#'
#' example("DESeq")
#' results(dds)
#' resultsNames(dds)
#' dds <- removeResults(dds)
#'
#' @rdname results
#' @aliases results resultsNames removeResults
#' @export
results <- function(object, name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  if (!"results" %in% mcols(mcols(object))$type) {
    stop("cannot find results columns in object, first call 'DESeq','nbinomWaldTest', or 'nbinomLRT'")
  }
  if (paste0("WaldPvalue_",name) %in% names(mcols(object))) {
    test <- "Wald"
  } else if ("LRTPvalue" %in% names(mcols(object))) {
    test <- "LRT"
  } else {
    stop("cannot find appropriate results, for available names call 'resultsNames(object)'")
  }
  log2FoldChange <- coef(object, name)
  lfcSE <- coefSE(object, name)
  pvalue <- pvalues(object,test,name)
  padj <- padj(object,test,name)
  res <- cbind(mcols(object)["baseMean"],
               log2FoldChange,
               lfcSE,
               pvalue,
               padj)
  names(res) <- c("baseMean","log2FoldChange","lfcSE","pvalue","padj")
  rownames(res) <- rownames(object)
  res
}

#' @rdname results
#' @export
resultsNames <- function(object) {
  names(mcols(object)[grep("log2 fold change",mcols(mcols(object))$description)])
}

#' @rdname results
#' @export
removeResults <- function(object) {
  resCols <- mcols(mcols(object))$type == "results"
  if (sum(resCols,na.rm=TRUE) > 0) {
    mcols(object) <- mcols(object)[,!resCols]
  }
  return(object)
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
  meanVarZero <- DataFrame(baseMean = rowMeans(counts(object,normalized=TRUE)),
                           baseVar = rowVars(counts(object,normalized=TRUE)),
                           allZero = rowSums(counts(object)) == 0)
  mcols(meanVarZero) <- DataFrame(type = rep("intermediate",ncol(meanVarZero)),
                                  description = c("the base mean over all rows",
                                    "the base variance over all rows",
                                    "all counts in a row are zero"))
  mcols(object) <- cbind(mcols(object),meanVarZero)
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
         warning( "dispersion fit did not converge" )
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
#
# return a list of results, with coefficients and standard
# errors on the log2 scale
fitNbinomGLMs <- function(object, modelMatrix, modelFormula, alpha_hat, lambda,
                          renameCols=TRUE, betaTol=1e-8, maxit=100, useOptim=TRUE) {
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
  
  if ("Intercept" %in% modelMatrixNames) {
    beta_mat <- matrix(0, ncol=ncol(modelMatrix), nrow=nrow(object))
    beta_mat[,which(modelMatrixNames == "Intercept")] <- log(mcols(object)$baseMean)
  } else {
    beta_mat <- matrix(1, ncol=ncol(modelMatrix), nrow=nrow(object))
  }
  
  # set a wide prior for all coefficients
  if (missing(lambda)) {
    lambda <- rep(1e-6, ncol(modelMatrix))
  }

  # here we convert from the log2 scale of the betas
  # to the log scale used in the C++ code
  # lambda is the inverse of the variance
  # so we use the squared inverse of the typical
  # conversion factor log(2)
  lambdaLogScale <- lambda * log2(exp(1))^2
  
  betaRes <- fitBeta(ySEXP = counts(object), xSEXP = modelMatrix,
                     nfSEXP = normalizationFactors,
                     alpha_hatSEXP = alpha_hat, beta_matSEXP = beta_mat,
                     lambdaSEXP = lambdaLogScale,
                     tolSEXP = betaTol, maxitSEXP = maxit)
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
  
  if (length(rowsForOptim) > 0) {
    x <- modelMatrix
    large <- 30
    for (row in rowsForOptim) {
      if (rowStable[row]) {
        betaRow <- betaMatrix[row,]
      } else {
        betaRow <- beta_mat[row,]
      }
      betaRow <- pmin(pmax(betaRow, -large), large)
      nf <- normalizationFactors[row,]
      k <- counts(object)[row,]
      alpha <- alpha_hat[row]
      objectiveFn <- function(p) {
        mu_row <- as.numeric(nf * 2^(x %*% p))
        prior <- sum(dnorm(p,0,sqrt(1/lambda),log=TRUE))
        logLike <- sum(dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE))
        -1 * (logLike + prior)
      }
      o <- optim(betaRow, objectiveFn, method="L-BFGS-B",
                 lower=-large, upper=large)
      ridge <- diag(lambda)
      # if we converged, change betaConv to TRUE
      if (o$convergence == 0) {
        betaConv[row] <- TRUE
      }
      # with or without convergence, store the estimate from optim
      betaMatrix[row,] <- o$par
      # calculate the standard errors
      mu_row <- as.numeric(nf * exp(x %*% (log(2)*o$par)))
      w <- diag((mu_row^-1 + alpha)^-1)
      xtwx <- t(x) %*% w %*% x
      xtwxRidgeInv <- solve(xtwx + ridge)
      betaSE[row,] <- log2(exp(1))*sqrt(diag(xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv))
      # store the new mu vector
      mu[row,] <- mu_row
      logLike[row] <- sum(dnbinom(k, mu_row, size=1/alpha, log=TRUE))
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
# beta_matSEXP n by k matrix of the initial estimates for the betas
# lambdaSEXP k length vector of the ridge values
# tolSEXP tolerance for convergence in estimates
# maxitSEXP maximum number of iterations
#
# return a list with elements: beta_mat, beta_var_mat, iter.
# Note: at this level the betas are on the natural log scale
fitBeta <- function (ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, beta_matSEXP, lambdaSEXP, tolSEXP, maxitSEXP) {
  # test for any NAs in arguments
  arg.names <- names(formals(fitBeta))
  na.test <- sapply(list(ySEXP, xSEXP, nfSEXP, alpha_hatSEXP, beta_matSEXP, lambdaSEXP, tolSEXP, maxitSEXP), function(x) any(is.na(x)))
  if (any(na.test)) {
    stop(paste("in call to fitBeta, the following arguments contain NA:",paste(arg.names[na.test],collapse=", ")))
  }
  .Call("fitBeta", ySEXP=ySEXP, xSEXP=xSEXP, nfSEXP=nfSEXP,
        alpha_hatSEXP=alpha_hatSEXP, beta_matSEXP=beta_matSEXP,
        lambdaSEXP=lambdaSEXP, tolSEXP=tolSEXP, maxitSEXP=maxitSEXP,
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
buildMatrixWithNARows <- function(m, NArows) {
  mFull <- matrix(NA, ncol=ncol(m), nrow=length(NArows))
  mFull[!NArows,] <- m
  mFull
}


# convenience function for breaking up matrices
# by column and preserving column names
matrixToList <- function(m) {
  l <- split(m, col(m))
  names(l) <- colnames(m)
  l
}

# convenience function to guess the name of the last coefficient
# in the model matrix, unless specified this will be used for
# plots and accessor functions
lastCoefName <- function(object) {
  resNms <- resultsNames(object)
  resNms[length(resNms)]
}


# functions to get coef, coefSE, pvalues and padj from mcols(object)
coef <- function(object,name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  mcols(object)[name]
}
coefSE <- function(object,name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  mcols(object)[paste0("SE_",name)]
}
pvalues <- function(object,test="Wald",name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  if (test == "Wald") {
    return(mcols(object)[paste0("WaldPvalue_",name)])
  } else if (test == "LRT") {
    return(mcols(object)["LRTPvalue"])
  } else {
    stop("unknown test")
  }
}
padj <- function(object,test="Wald",name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  if (test == "Wald") {
    return(mcols(object)[paste0("WaldAdjPvalue_",name)])
  } else if (test == "LRT") {
    return(mcols(object)["LRTAdjPvalue"])
  } else {
    stop("unknown test")
  }
}

# convenience function to make more descriptive names
# for factor variables
renameModelMatrixColumns <- function(modelMatrixNames, data, design) {
  designVars <- all.vars(formula(design))
  designVarsClass <- sapply(designVars, function(v) class(data[[v]]))
  factorVars <- designVars[designVarsClass == "factor"]
  colNamesFrom <- do.call(c,lapply(factorVars, function(v) paste0(v,levels(data[[v]])[-1])))
  colNamesTo <- do.call(c,lapply(factorVars, function(v) paste0(v,"_",levels(data[[v]])[-1],"_vs_",levels(data[[v]])[1])))
  data.frame(from=colNamesFrom,to=colNamesTo,stringsAsFactors=FALSE)
}

calculateCooksDistance <- function(object, H, p) {
  V <- assays(object)[["mu"]] + mcols(object)$dispFit * assays(object)[["mu"]]^2
  PearsonResSq <- (counts(object) - assays(object)[["mu"]])^2 / V
  cooks <- PearsonResSq / p  * H / (1 - H)^2
  cooks
}

leaveOneOutFullRank <- function(modelMatrix) {
  sapply(seq_len(nrow(modelMatrix)), function(i) qr(modelMatrix[-i,])$rank) == ncol(modelMatrix)
}
