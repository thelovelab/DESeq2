#' Extract results from a DESeq analysis
#'
#' \code{results} extracts a result table from a DESeq analysis giving base means across samples,
#' log2 fold changes, standard errors, test statistics, p-values and adjusted p-values;
#' \code{resultsNames} returns the names of the estimated effects (coefficents) of the model;
#' \code{removeResults} returns a \code{DESeqDataSet} object with results columns removed.
#'
#' The results table when printed will provide the information about
#' the comparison, e.g. "log2 fold change (MAP): condition treated vs untreated", meaning
#' that the estimates are of log2(treated / untreated), as would be returned by
#' \code{contrast=c("condition","treated","untreated")}.
#' Multiple results can be returned for analyses beyond a simple two group comparison,
#' so \code{results} takes arguments \code{contrast} and \code{name} to help
#' the user pick out the comparisons of interest for printing a results table.
#' The use of the \code{contrast} argument is recommended for exact specification
#' of the levels which should be compared and their order.
#' 
#' If \code{results} is run without specifying \code{contrast} or \code{name},
#' it will return the comparison of the last level of the last variable in the
#' design formula over the first level of this variable. For example, for a simple two-group
#' comparison, this would return the log2 fold changes of the second group over the
#' first group (the reference level). Please see examples below and in the vignette. 
#'
#' The argument \code{contrast} can be used to generate results tables for
#' any comparison of interest, for example, the log2 fold change between
#' two levels of a factor, and its usage is described below. It can also
#' accomodate more complicated numeric comparisons.
#' Note that \code{contrast} will set to 0 the estimated LFC in a
#' comparison of two groups, where all of the counts in the two groups
#' are equal to 0 (while other groups have positive counts), while
#' \code{name} will not automatically set these LFC to 0.
#' The test statistic used for a contrast is:
#'
#' \deqn{ c^t \beta / \sqrt{c^t \Sigma c } }{ c' beta / sqrt( c' Sigma c ) }
#'
#' The argument \code{name} can be used to generate results tables for
#' individual effects, which must be individual elements of \code{resultsNames(object)}.
#' These individual effects could represent continuous covariates, effects
#' for individual levels, or individual interaction effects.
#' 
#' Information on the comparison which was used to build the results table,
#' and the statistical test which was used for p-values (Wald test or likelihood ratio test)
#' is stored within the object returned by \code{results}. This information is in
#' the metadata columns of the results table, which is accessible by calling \code{mcols}
#' on the \code{\link{DESeqResults}} object returned by \code{results}.
#'
#' On p-values:
#' 
#' By default, independent filtering is performed to select a set of genes
#' for multiple test correction which maximizes the number of adjusted
#' p-values less than a given critical value \code{alpha} (by default 0.1).
#' See the reference in this man page for details on independent filtering.
#' The filter used for maximizing the number of rejections is the mean
#' of normalized counts for all samples in the dataset.
#' Several arguments from the \code{filtered_p} function of
#' the genefilter package (used within the \code{results} function)
#' are provided here to control the independent filtering behavior.
#' (Note \code{filtered_p} R code is now copied into DESeq2
#' package to avoid gfortran requirements.)
#' In DESeq2 version >= 1.10, the threshold that is chosen is
#' the lowest quantile of the filter for which the
#' number of rejections is close to the peak of a curve fit
#' to the number of rejections over the filter quantiles.
#' 'Close to' is defined as within 1 residual standard deviation.
#' The adjusted p-values for the genes which do not pass the filter threshold
#' are set to \code{NA}. 
#' 
#' By default, \code{results} assigns a p-value of \code{NA}
#' to genes containing count outliers, as identified using Cook's distance.
#' See the \code{cooksCutoff} argument for control of this behavior.
#' Cook's distances for each sample are accessible as a matrix "cooks"
#' stored in the \code{assays()} list. This measure is useful for identifying rows where the
#' observed counts might not fit to a Negative Binomial distribution.
#'
#' For analyses using the likelihood ratio test (using \code{\link{nbinomLRT}}),
#' the p-values are determined solely by the difference in deviance between
#' the full and reduced model formula. A single log2 fold change is printed
#' in the results table for consistency with other results table outputs,
#' however the test statistic and p-values may nevertheless involve
#' the testing of one or more log2 fold changes.
#' Which log2 fold change is printed in the results table can be controlled
#' using the \code{name} argument, or by default this will be the estimated
#' coefficient for the last element of \code{resultsNames(object)}.
#'
#' If \code{useT=TRUE} was specified when running \code{DESeq} or \code{nbinomWaldTest},
#' then the p-value generated by \code{results} will also make use of the
#' t distribution for the Wald statistic, using the degrees of freedom
#' in \code{mcols(object)$tDegreesFreedom}.
#' 
#' @references Richard Bourgon, Robert Gentleman, Wolfgang Huber: Independent
#' filtering increases detection power for high-throughput experiments.
#' PNAS (2010), \url{http://dx.doi.org/10.1073/pnas.0914005107}
#' 
#' @param object a DESeqDataSet, on which one
#' of the following functions has already been called:
#' \code{\link{DESeq}}, \code{\link{nbinomWaldTest}}, or \code{\link{nbinomLRT}}
#' @param contrast this argument specifies what comparison to extract from
#' the \code{object} to build a results table. one of either:
#' \itemize{
#'  \item a character vector with exactly three elements:
#' the name of a factor in the design formula,
#' the name of the numerator level for the fold change,
#' and the name of the denominator level for the fold change
#' (simplest case)
#'  \item a list of 2 character vectors: the names of the fold changes
#' for the numerator, and the names of the fold changes
#' for the denominator.
#' these names should be elements of \code{resultsNames(object)}.
#' if the list is length 1, a second element is added which is the
#' empty character vector, \code{character()}.
#' (more general case, can be to combine interaction terms and main effects)
#'  \item a numeric contrast vector with one element
#' for each element in \code{resultsNames(object)} (most general case)
#' }
#' If specified, the \code{name} argument is ignored.
#' @param name the name of the individual effect (coefficient) for
#' building a results table. Use this argument rather than \code{contrast}
#' for continuous variables, individual effects or for individual interaction terms.
#' The value provided to \code{name} must be an element of \code{resultsNames(object)}.
#' @param lfcThreshold a non-negative value which specifies a log2 fold change
#' threshold. The default value is 0, corresponding to a test that
#' the log2 fold changes are equal to zero. The user can
#' specify the alternative hypothesis using the \code{altHypothesis} argument,
#' which defaults to testing
#' for log2 fold changes greater in absolute value than a given threshold.
#' If \code{lfcThreshold} is specified,
#' the results are for Wald tests, and LRT p-values will be overwritten.
#' @param altHypothesis character which specifies the alternative hypothesis,
#' i.e. those values of log2 fold change which the user is interested in
#' finding. The complement of this set of values is the null hypothesis which
#' will be tested. If the log2 fold change specified by \code{name}
#' or by \code{contrast} is written as \eqn{ \beta }{ beta }, then the possible values for
#' \code{altHypothesis} represent the following alternate hypotheses:
#' \itemize{
#' \item greaterAbs: \eqn{|\beta| > \textrm{lfcThreshold} }{ |beta| > lfcThreshold },
#' and p-values are two-tailed
#' \item lessAbs: \eqn{ |\beta| < \textrm{lfcThreshold} }{ |beta| < lfcThreshold },
#' p-values are the maximum of the upper and lower tests.
#' The Wald statistic given is positive, an SE-scaled distance from the closest boundary
#' \item greater: \eqn{ \beta > \textrm{lfcThreshold} }{ beta > lfcThreshold }
#' \item less: \eqn{ \beta < -\textrm{lfcThreshold} }{ beta < -lfcThreshold }
#' }
#' @param listValues only used if a list is provided to \code{contrast}:
#' a numeric of length two: the log2 fold changes in the list are multiplied by these values.
#' the first number should be positive and the second negative. 
#' by default this is \code{c(1,-1)}
#' @param cooksCutoff theshold on Cook's distance, such that if one or more
#' samples for a row have a distance higher, the p-value for the row is
#' set to NA. The default cutoff is the .99 quantile of the F(p, m-p) distribution,
#' where p is the number of coefficients being fitted and m is the number of samples.
#' Set to \code{Inf} or \code{FALSE} to disable the resetting of p-values to NA.
#' Note: this test excludes the Cook's distance of samples belonging to experimental
#' groups with only 2 samples.
#' @param independentFiltering logical, whether independent filtering should be
#' applied automatically
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering (by default 0.1). If the adjusted p-value cutoff (FDR) will be a
#' value other than 0.1, \code{alpha} should be set to that value.
#' @param filter the vector of filter statistics over which the independent
#' filtering will be optimized. By default the mean of normalized counts is used.
#' @param theta the quantiles at which to assess the number of rejections
#' from independent filtering
#' @param pAdjustMethod the method to use for adjusting p-values, see \code{?p.adjust}
#' @param filterFun an optional custom function for performing independent filtering
#' and p-value adjustment, with arguments \code{res} (a DESeqResults object),
#' \code{filter} (the quantitity for filtering tests),
#' \code{alpha} (the target FDR),
#' \code{pAdjustMethod}. This function should return a DESeqResults object
#' with a \code{padj} column.
#' @param format character, either \code{"DataFrame"},
#' \code{"GRanges"}, or \code{"GRangesList"},
#' whether the results should be printed as a \code{\link{DESeqResults}} DataFrame,
#' or if the results DataFrame should be attached as metadata columns to
#' the \code{GRanges} or \code{GRangesList} \code{rowRanges} of the \code{DESeqDataSet}.
#' If the \code{rowRanges} is a \code{GRangesList}, and \code{GRanges} is requested, 
#' the range of each gene will be returned
#' @param saveCols character or numeric vector, the columns of
#' \code{mcols(object)} to pass into the \code{results} output
#' @param test this is automatically detected internally if not provided.
#' the one exception is after \code{nbinomLRT} has been run, \code{test="Wald"}
#' will generate Wald statistics and Wald test p-values.
#' @param addMLE if \code{betaPrior=TRUE} was used (non-default),
#' this logical argument specifies if the "unshrunken" maximum likelihood estimates (MLE)
#' of log2 fold change should be added as a column to the results table (default is FALSE).
#' This argument is preserved for backward compatability, as now \code{betaPrior=TRUE}
#' by default and the recommended pipeline is
#' to generate shrunken MAP estimates using \code{\link{lfcShrink}}.
#' This argument functionality is only implemented for \code{contrast}
#' specified as three element character vectors.
#' @param tidy whether to output the results table with rownames as a first column 'row'.
#' the table will also be coerced to \code{data.frame}
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see next argument \code{BPPARAM}
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}.
#' If not specified, the parameters last registered with
#' \code{\link{register}} will be used.
#' @param minmu lower bound on the estimated count (used when calculating contrasts)
#' 
#' @return For \code{results}: a \code{\link{DESeqResults}} object, which is
#' a simple subclass of DataFrame. This object contains the results columns:
#' \code{baseMean}, \code{log2FoldChange}, \code{lfcSE}, \code{stat},
#' \code{pvalue} and \code{padj},
#' and also includes metadata columns of variable information.
#' The \code{lfcSE} gives the standard error of the \code{log2FoldChange}.
#' For the Wald test, \code{stat} is the Wald statistic: the \code{log2FoldChange}
#' divided by \code{lfcSE}, which is compared to a standard Normal distribution
#' to generate a two-tailed \code{pvalue}. For the likelihood ratio test (LRT),
#' \code{stat} is the difference in deviance between the reduced model and the full model,
#' which is compared to a chi-squared distribution to generate a \code{pvalue}.
#'
#' For \code{resultsNames}: the names of the columns available as results,
#' usually a combination of the variable name and a level
#'
#' For \code{removeResults}: the original \code{DESeqDataSet} with results metadata columns removed
#'
#' @seealso \code{\link{DESeq}}, \code{\link{lfcShrink}}
#'
#' @examples
#'
#' ## Example 1: two-group comparison
#' 
#' dds <- makeExampleDESeqDataSet(m=4)
#' 
#' dds <- DESeq(dds)
#' res <- results(dds, contrast=c("condition","B","A"))
#'
#' # with more than two groups, the call would look similar, e.g.:
#' # results(dds, contrast=c("condition","C","A"))
#' # etc.
#' 
#' ## Example 2: two conditions, two genotypes, with an interaction term
#' 
#' dds <- makeExampleDESeqDataSet(n=100,m=12)
#' dds$genotype <- factor(rep(rep(c("I","II"),each=3),2))
#' 
#' design(dds) <- ~ genotype + condition + genotype:condition
#' dds <- DESeq(dds) 
#' resultsNames(dds)
#'
#' # the condition effect for genotype I (the main effect)
#' results(dds, contrast=c("condition","B","A"))
#'
#' # the condition effect for genotype II
#' # this is, by definition, the main effect *plus* the interaction term
#' # (the extra condition effect in genotype II compared to genotype I).
#' results(dds, list( c("condition_B_vs_A","genotypeII.conditionB") ))
#' 
#' # the interaction term, answering: is the condition effect *different* across genotypes?
#' results(dds, name="genotypeII.conditionB")
#'  
#' ## Example 3: two conditions, three genotypes
#'
#' # ~~~ Using interaction terms ~~~
#' 
#' dds <- makeExampleDESeqDataSet(n=100,m=18)
#' dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
#' design(dds) <- ~ genotype + condition + genotype:condition
#' dds <- DESeq(dds)
#' resultsNames(dds)
#'
#' # the condition effect for genotype I (the main effect)
#' results(dds, contrast=c("condition","B","A"))
#'
#' # the condition effect for genotype III.
#' # this is the main effect *plus* the interaction term
#' # (the extra condition effect in genotype III compared to genotype I).
#' results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))
#'  
#' # the interaction term for condition effect in genotype III vs genotype I.
#' # this tests if the condition effect is different in III compared to I
#' results(dds, name="genotypeIII.conditionB")
#' 
#' # the interaction term for condition effect in genotype III vs genotype II.
#' # this tests if the condition effect is different in III compared to II
#' results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))
#'
#' # Note that a likelihood ratio could be used to test if there are any
#' # differences in the condition effect between the three genotypes.
#' 
#' # ~~~ Using a grouping variable ~~~
#' 
#' # This is a useful construction when users just want to compare
#' # specific groups which are combinations of variables.
#' 
#' dds$group <- factor(paste0(dds$genotype, dds$condition))
#' design(dds) <- ~ group
#' dds <- DESeq(dds)
#' resultsNames(dds)
#'
#' # the condition effect for genotypeIII
#' results(dds, contrast=c("group", "IIIB", "IIIA"))
#' 
#' @rdname results
#' @aliases results resultsNames removeResults
#' @export
results <- function(object, contrast, name, 
                    lfcThreshold=0,
                    altHypothesis=c("greaterAbs","lessAbs","greater","less"),
                    listValues=c(1,-1),
                    cooksCutoff,
                    independentFiltering=TRUE,
                    alpha=0.1, filter, theta,
                    pAdjustMethod="BH",
                    filterFun,
                    format=c("DataFrame","GRanges","GRangesList"),
                    saveCols=NULL,
                    test,
                    addMLE=FALSE,
                    tidy=FALSE,
                    parallel=FALSE, BPPARAM=bpparam(), 
                    minmu=0.5) {

  stopifnot(is(object, "DESeqDataSet"))
  
  # match args
  format <- match.arg(format, choices=c("DataFrame", "GRanges","GRangesList"))
  altHypothesis <- match.arg(altHypothesis, choices=c("greaterAbs","lessAbs","greater","less"))
  if (!missing(test)) {
    test <- match.arg(test, choices=c("Wald","LRT"))
  }
  
  ### initial argument testing ###
  
  stopifnot(lfcThreshold >= 0)
  stopifnot(length(lfcThreshold)==1)
  stopifnot(length(alpha)==1)
  stopifnot(alpha > 0 & alpha < 1)
  stopifnot(length(pAdjustMethod)==1)
  stopifnot(length(listValues)==2 & is.numeric(listValues))
  stopifnot(listValues[1] > 0 & listValues[2] < 0)
  if (!"results" %in% mcols(mcols(object))$type) {
    stop("couldn't find results. you should first run DESeq()")
  }
  if (missing(test)) {
    test <- attr(object, "test")
  } else if (test == "Wald" & attr(object, "test") == "LRT") {
    # initially test was LRT, now need to add Wald statistics and p-values
    object <- makeWaldTest(object)
  } else if (test == "LRT" & attr(object, "test") == "Wald") {
    stop("the LRT requires the user run nbinomLRT or DESeq(dds,test='LRT')")
  }
  if (lfcThreshold == 0 & altHypothesis == "lessAbs") {
    stop("when testing altHypothesis='lessAbs', set the argument lfcThreshold to a positive value")
  }
  
  if (addMLE) {
    if (!attr(object,"betaPrior")) {
      stop("addMLE=TRUE is only for when a beta prior was used.
otherwise, the log2 fold changes are already MLE")
    }
    if (!missing(name) & missing(contrast)) {
      stop("addMLE=TRUE should be used by providing character vector
of length 3 to 'contrast' instead of using 'name'")
    }
  }
  
  if (format == "GRanges" & is(rowRanges(object),"GRangesList")) {
    if (any(elementNROWS(rowRanges(object)) == 0)) {
      stop("rowRanges is GRangesList and one or more GRanges have length 0. Use format='DataFrame' or 'GRangesList'")
    }
  }

  if (!is.null(saveCols)) {
    if (is(saveCols, "character"))
      stopifnot(all(saveCols %in% colnames(mcols(object))))
    if (is(saveCols, "numeric")) {
      stopifnot(saveCols == round(saveCols))
      stopifnot(min(saveCols) > 0)
      stopifnot(max(saveCols) <= ncol(mcols(object)))
    }
  }

  if (!missing(contrast)) {
    if (attr(object,"modelMatrixType") == "user-supplied" & is.character(contrast)) {
      stop("only list- and numeric-type contrasts are supported for user-supplied model matrices")
    }
  }

  if (is(design(object), "formula")) {
    hasIntercept <- attr(terms(design(object)),"intercept") == 1
    isExpanded <- attr(object, "modelMatrixType") == "expanded"
    termsOrder <- attr(terms.formula(design(object)),"order")
    # if no intercept was used or an expanded model matrix was used, 
    # and neither 'contrast' nor 'name' were specified,
    # and no interactions...
    # then we create the result table: last / first level for last variable
    if ((test == "Wald") & (isExpanded | !hasIntercept) & missing(contrast) & missing(name) & all(termsOrder < 2)) {
      designVars <- all.vars(design(object))
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      if (is.factor(lastVar)) {
        nlvls <- nlevels(lastVar)
        contrast <- c(lastVarName, levels(lastVar)[nlvls], levels(lastVar)[1])
      }
    }
  }
  
  if (missing(name)) {
    name <- lastCoefName(object)
  } else { 
    if (length(name) != 1 | !is.character(name)) {
      stop("the argument 'name' should be a character vector of length 1")
    }
  }

  ### done with input argument testing ###
  
  WaldResults <- paste0("WaldPvalue_",name) %in% names(mcols(object))
  LRTResults <- "LRTPvalue" %in% names(mcols(object))

  # this will be used in cleanContrast, and in the lfcThreshold chunks below
  useT <- "tDegreesFreedom" %in% names(mcols(object))
  
  # if performing a contrast call the function cleanContrast()
  if (!missing(contrast)) {
    resNames <- resultsNames(object)
    # do some arg checking/cleaning
    contrast <- checkContrast(contrast, resNames)

    ### cleanContrast call ###   
    # need to go back to C++ code in order to build the beta covariance matrix
    # then this is multiplied by the numeric contrast to get the Wald statistic.
    # with 100s of samples, this can get slow, so offer parallelization
    if (!parallel) {
      res <- cleanContrast(object, contrast, expanded=isExpanded, listValues=listValues,
                           test=test, useT=useT, minmu=minmu)
    } else if (parallel) {
      # parallel execution
      nworkers <- getNworkers(BPPARAM)
      idx <- factor(sort(rep(seq_len(nworkers),length.out=nrow(object))))
      res <- do.call(rbind, bplapply(levels(idx), function(l) {
        cleanContrast(object[idx == l,,drop=FALSE], contrast,
                      expanded=isExpanded, listValues=listValues,
                      test=test, useT=useT, minmu=minmu)
      }, BPPARAM=BPPARAM))
    }

  } else {
    # if not performing a contrast
    # pull relevant columns from mcols(object)
    log2FoldChange <- getCoef(object, name)
    lfcSE <- getCoefSE(object, name)
    stat <- getStat(object, test, name)
    pvalue <- getPvalue(object, test, name)
    res <- cbind(mcols(object)["baseMean"],log2FoldChange,lfcSE,stat,pvalue)
    names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
  }
  
  rownames(res) <- rownames(object)

  # add unshrunken MLE coefficients to the results table
  if (addMLE) {
    if (is.numeric(contrast)) stop("addMLE only implemented for: contrast=c('condition','B','A')")
    if (is.list(contrast)) stop("addMLE only implemented for: contrast=c('condition','B','A')")
    res <- cbind(res, mleContrast(object, contrast))
    res <- res[,c("baseMean","log2FoldChange","lfcMLE","lfcSE","stat","pvalue")]
    # if an all zero contrast, also zero out the lfcMLE
    res$lfcMLE[ which(res$log2FoldChange == 0 & res$stat == 0) ] <- 0
  }
  
  # only if we need to generate new p-values
  if ( !(lfcThreshold == 0 & altHypothesis == "greaterAbs") ) {
    if (test == "LRT") {
      stop("tests of log fold change above or below a theshold must be Wald tests.")
    }
    # check requirement if betaPrior was set to FALSE
    if (altHypothesis == "lessAbs" & attr(object, "betaPrior")) {
      stop("testing altHypothesis='lessAbs' requires setting the DESeq() argument betaPrior=FALSE")
    }
    # easier to read
    LFC <- res$log2FoldChange
    SE <- res$lfcSE
    T <- lfcThreshold

    if (useT) {
      df <- mcols(object)$tDegreesFreedom
      pfunc <- function(q) pt(q, df=df, lower.tail=FALSE)
    } else {
      pfunc <- function(q) pnorm(q, lower.tail=FALSE)
    }
    
    if (altHypothesis == "greaterAbs") {
      newStat <- sign(LFC) * pmax((abs(LFC) - T)/SE, 0)
      newPvalue <- pmin(1, 2 * pfunc((abs(LFC) - T)/SE))
    } else if (altHypothesis == "lessAbs") {
      newStatAbove <- pmax((T - LFC)/SE, 0)
      pvalueAbove <- pfunc((T - LFC)/SE)
      newStatBelow <- pmax((LFC + T)/SE, 0)
      pvalueBelow <- pfunc((LFC + T)/SE)
      newStat <- pmin(newStatAbove, newStatBelow)
      newPvalue <- pmax(pvalueAbove, pvalueBelow)
    } else if (altHypothesis == "greater") {
      newStat <- pmax((LFC - T)/SE, 0)
      newPvalue <- pfunc((LFC - T)/SE)
    } else if (altHypothesis == "less") {
      newStat <- pmin((LFC + T)/SE, 0)
      newPvalue <- pfunc((-T - LFC)/SE)
    }
    res$stat <- newStat
    res$pvalue <- newPvalue
  }
  
  # calculate Cook's cutoff
  m <- nrow(attr(object,"dispModelMatrix"))
  p <- ncol(attr(object,"dispModelMatrix"))
  
  defaultCutoff <- qf(.99, p, m - p)
  if (missing(cooksCutoff)) {
    cooksCutoff <- defaultCutoff
  }
  stopifnot(length(cooksCutoff)==1)
  if (is.logical(cooksCutoff) & cooksCutoff) {
    cooksCutoff <- defaultCutoff
  }
  
  # apply cutoff based on maximum Cook's distance
  performCooksCutoff <- (is.numeric(cooksCutoff) | cooksCutoff) 
  if (performCooksCutoff) {
    cooksOutlier <- mcols(object)$maxCooks > cooksCutoff
    
    ### BEGIN heuristic to avoid filtering genes with low count outliers
    # as according to Cook's cutoff. only for two group designs.
    # dont filter if three or more counts are larger
    if (any(cooksOutlier,na.rm=TRUE) & is(design(object), "formula")) {
      designVars <- all.vars(design(object))
      if (length(designVars) == 1) {
        var <- colData(object)[[designVars]]
        if (is(var, "factor") && nlevels(var) == 2) {
          dontFilter <- logical(sum(cooksOutlier,na.rm=TRUE))
          for (i in seq_along(dontFilter)) {
            # index along rows of object
            ii <- which(cooksOutlier)[i] 
            # count for the outlier with max cooks
            outCount <- counts(object)[ii,which.max(assays(object)[["cooks"]][ii,])]
            # if three or more counts larger than the outlier
            if (sum(counts(object)[ii,] > outCount) >= 3) {
              # don't filter out the p-value for that gene
              dontFilter[i] <- TRUE
            }
          }
          # reset the outlier status for these genes
          cooksOutlier[which(cooksOutlier)][dontFilter] <- FALSE
        }
      }
    } ### END heuristic ###

    res$pvalue[which(cooksOutlier)] <- NA
  }

  # if original baseMean was positive, but now zero due to replaced counts, fill in results
  if ( sum(mcols(object)$replace, na.rm=TRUE) > 0) {
    nowZero <- which(mcols(object)$replace & mcols(object)$baseMean == 0)
    res$log2FoldChange[nowZero] <- 0
    if (addMLE) { res$lfcMLE[nowZero] <- 0 }
    res$lfcSE[nowZero] <- 0
    res$stat[nowZero] <- 0
    res$pvalue[nowZero] <- 1
  }

  # add prior information
  deseq2.version <- packageVersion("DESeq2")
  if (!attr(object,"betaPrior")) {
    priorInfo <- list(type="none", package="DESeq2", version=deseq2.version)
  } else {
    betaPriorVar <- attr(object, "betaPriorVar")
    priorInfo <- list(type="normal", package="DESeq2", version=deseq2.version,
                      betaPriorVar=betaPriorVar)
  }
  
  # make results object
  deseqRes <- DESeqResults(res, priorInfo=priorInfo)
  
  # p-value adjustment
  if (missing(filterFun)) {
    deseqRes <- pvalueAdjustment(deseqRes, independentFiltering, filter,
                                 theta, alpha, pAdjustMethod)
  } else {
    deseqRes <- filterFun(deseqRes, filter, alpha, pAdjustMethod)
  }

  # stash lfcThreshold
  metadata(deseqRes)[["lfcThreshold"]] <- lfcThreshold
  
  # remove rownames and attach as a new column, 'row'
  if (tidy) {
    colnms <- colnames(deseqRes)
    deseqRes$row <- rownames(deseqRes)
    mcols(deseqRes,use.names=TRUE)["row","type"] <- "results"
    mcols(deseqRes,use.names=TRUE)["row","description"] <- "row names"
    deseqRes <- deseqRes[,c("row",colnms)]
    rownames(deseqRes) <- NULL
    deseqRes <- as.data.frame(deseqRes)
  }

  # return DataFrame, GRanges or GRangesList
  out <- resultsFormatSwitch(object=object, res=deseqRes, format=format, saveCols=saveCols)
  return(out)
}

#' @rdname results
#' @export
resultsNames <- function(object) {
  names(mcols(object))[grep("log2 fold change",mcols(mcols(object))$description)]
}

#' @rdname results
#' @export
removeResults <- function(object) {
  resCols <- mcols(mcols(object))$type == "results"
  if (sum(resCols,na.rm=TRUE) > 0) {
    mcols(object) <- mcols(object)[,-which(resCols),drop=FALSE]
  }
  return(object)
}


###########################################################
# unexported functons 
###########################################################

pvalueAdjustment <- function(res, independentFiltering, filter,
                             theta, alpha, pAdjustMethod) {
  # perform independent filtering
  if (independentFiltering) {
    if (missing(filter)) {
      filter <- res$baseMean
    }
    if (missing(theta)) {
      lowerQuantile <- mean(filter == 0)
      if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
      theta <- seq(lowerQuantile, upperQuantile, length=50)
    }

    # do filtering using genefilter
    stopifnot(length(theta) > 1)
    stopifnot(length(filter) == nrow(res))
    filtPadj <- filtered_p(filter=filter, test=res$pvalue,
                           theta=theta, method=pAdjustMethod) 
    numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
    # prevent over-aggressive filtering when all genes are null,
    # by requiring the max number of rejections is above a fitted curve.
    # If the max number of rejection is not greater than 10, then don't
    # perform independent filtering at all.
    lo.fit <- lowess(numRej ~ theta, f=1/5)
    if (max(numRej) <= 10) {
      j <- 1
    } else { 
      residual <- if (all(numRej==0)) {
        0
      } else {
        numRej[numRej > 0] - lo.fit$y[numRej > 0]
      }
      thresh <- max(lo.fit$y) - sqrt(mean(residual^2))
      j <- if (any(numRej > thresh)) {
        which(numRej > thresh)[1]
      } else {
        1  
      }
    }
    # j <- which.max(numRej) # old method
    padj <- filtPadj[, j, drop=TRUE]
    cutoffs <- quantile(filter, theta)
    filterThreshold <- cutoffs[j]
    filterNumRej <- data.frame(theta=theta, numRej=numRej)
    filterTheta <- theta[j]

    metadata(res)[["filterThreshold"]] <- filterThreshold
    metadata(res)[["filterTheta"]] <- filterTheta
    metadata(res)[["filterNumRej"]] <- filterNumRej
    metadata(res)[["lo.fit"]] <- lo.fit
    metadata(res)[["alpha"]] <- alpha
  } else {
    # regular p-value adjustment
    # does not include those rows which were removed
    # by maximum Cook's distance
    padj <- p.adjust(res$pvalue,method=pAdjustMethod)  
  }
  res$padj <- padj
  
  # add metadata to padj column
  mcols(res)$type[names(res)=="padj"] <- "results"
  mcols(res)$description[names(res)=="padj"] <- paste(pAdjustMethod,"adjusted p-values")
  
  res
}

# function copied from `genefilter` package to avoid gfortran requirement
filtered_p <- function( filter, test, theta, data, method = "none" ) {
  if ( is.function( filter ) )
    U1 <- filter( data )
  else
    U1 <- filter
  cutoffs <- quantile( U1, theta )
  result <- matrix( NA_real_, length( U1 ), length( cutoffs ) )
  colnames( result ) <- names( cutoffs )
  for ( i in 1:length( cutoffs ) ) {    
    use <- U1 >= cutoffs[i]
    if( any( use ) ) {
      if( is.function( test ) )
        U2 <- test( data[use,] )
      else
        U2 <- test[use]
      result[use,i] <- p.adjust( U2, method )
    }
  }
  return( result )
}

# two low-level functions used by results() to perform contrasts
#
# getContrast takes a DESeqDataSet object
# and a numeric vector specifying a contrast
# and returns a vector of Wald statistics
# corresponding to the contrast.
#
# cleanContrast checks for the validity of
# the specified contrast (numeric or character vector)
# and turns character vector contrast into the appropriate
# numeric vector contrast
#
# results() calls cleanContrast() which calls getContrast()
#
# the formula used is:
# c' beta / sqrt( c' sigma c)
# where beta is the coefficient vector
# and sigma is the covariance matrix for beta
getContrast <- function(object, contrast, useT=FALSE, minmu) {
  if (missing(contrast)) {
    stop("must provide a contrast")
  }
  if (is.null(attr(object,"modelMatrix"))) {
    stop("was expecting a model matrix stored as an attribute of the DESeqDataSet")
  }
  modelMatrix <- attr(object, "modelMatrix")
  
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]
  normalizationFactors <- getSizeOrNormFactors(objectNZ)
  alpha_hat <- dispersions(objectNZ)
  coefColumns <- names(mcols(objectNZ))[grep("log2 fold change",mcols(mcols(object))$description)]
  # convert betas to log scale
  beta_mat <- log(2) * as.matrix(mcols(objectNZ)[,coefColumns,drop=FALSE])
  # convert beta prior variance to log scale
  lambda = 1/(log(2)^2 * attr(object,"betaPriorVar"))

  # check if DESeq() replaced outliers
  countsMatrix <- if ("replaceCounts" %in% assayNames(object)) {
    assays(objectNZ)[["replaceCounts"]]
  } else {
    counts(objectNZ)
  }

  # use weights if they are present in assays(object)
  if ("weights" %in% assayNames(object)) {
    useWeights <- TRUE
    weights <- assays(object)[["weights"]]
    stopifnot(all(weights >= 0))
    weights <- weights / apply(weights, 1, max)
  } else {
    useWeights <- FALSE
    weights <- matrix(1, nrow=nrow(object), ncol=ncol(object))
  }
  
  betaRes <- fitBeta(ySEXP = countsMatrix, xSEXP = modelMatrix,
                     nfSEXP = normalizationFactors,
                     alpha_hatSEXP = alpha_hat,
                     contrastSEXP = contrast,
                     beta_matSEXP = beta_mat,
                     lambdaSEXP = lambda,
                     weightsSEXP = weights,
                     useWeightsSEXP = useWeights,
                     tolSEXP = 1e-8, maxitSEXP = 0,
                     useQRSEXP=FALSE, # QR not relevant, fitting loop isn't entered
                     minmuSEXP=minmu)
  # convert back to log2 scale
  contrastEstimate <- log2(exp(1)) * betaRes$contrast_num
  contrastSE <- log2(exp(1)) * betaRes$contrast_denom
  contrastStatistic <- contrastEstimate / contrastSE
  if (useT) {
    stopifnot("tDegreesFreedom" %in% names(mcols(object)))
    df <- mcols(objectNZ)$tDegreesFreedom
    contrastPvalue <- 2*pt(abs(contrastStatistic),df=df,lower.tail=FALSE)
  } else {
    contrastPvalue <- 2*pnorm(abs(contrastStatistic),lower.tail=FALSE)
  }
  contrastList <- list(log2FoldChange=contrastEstimate,
                       lfcSE=contrastSE,
                       stat=contrastStatistic,
                       pvalue=contrastPvalue)
  contrastResults <- buildDataFrameWithNARows(contrastList,
                                              mcols(object)$allZero)
  names(contrastResults) <- c("log2FoldChange","lfcSE","stat","pvalue")
  contrastResults
}

# this function takes a desired contrast as specified by results(),
# performs checks, and then either returns the already existing contrast
# or generates the contrast by calling getContrast() using a numeric vector
cleanContrast <- function(object, contrast, expanded=FALSE, listValues, test, useT, minmu) {
  # get the names of columns in the beta matrix
  resNames <- resultsNames(object)
  # if possible, return pre-computed columns, which are
  # already stored in mcols(dds). this will be the case using
  # results() with 'name', or if expanded model matrices were not
  # run and the contrast contains the reference level as numerator or denominator

  resReady <- FALSE
  
  if (is.character(contrast)) {
    contrastFactor <- contrast[1]
    if (!contrastFactor %in% names(colData(object))) {
      stop(paste(contrastFactor,"should be the name of a factor in the colData of the DESeqDataSet"))
    }
    if (!is(colData(object)[[contrastFactor]], "factor")) {
      stop(paste(contrastFactor,"is not a factor, see ?results"))
    }
    contrastNumLevel <- contrast[2]
    contrastDenomLevel <- contrast[3]
    contrastBaseLevel <- levels(colData(object)[,contrastFactor])[1]
    
    # check for intercept
    hasIntercept <- attr(terms(design(object)),"intercept") == 1
    firstVar <- contrastFactor == all.vars(design(object))[1]

    # tricky case: if the design has no intercept, the factor is
    # not the first variable in the design, and one of the numerator or denominator
    # is the reference level, then the desired contrast is simply a coefficient (or -1 times)
    noInterceptPullCoef <- !hasIntercept & !firstVar &
      (contrastBaseLevel %in% c(contrastNumLevel, contrastDenomLevel))
    
    # case 1: standard model matrices: pull coef or build the appropriate contrast
    # coefficients names are of the form  "factor_level_vs_baselevel"
    # output: contrastNumColumn and contrastDenomColumn
    if (!expanded & (hasIntercept | noInterceptPullCoef)) {
      # use make.names() so the column names are
      # the same as created by DataFrame in mcols(object).
      contrastNumColumn <- make.names(paste0(contrastFactor,"_",contrastNumLevel,"_vs_",contrastBaseLevel))
      contrastDenomColumn <- make.names(paste0(contrastFactor,"_",contrastDenomLevel,"_vs_",contrastBaseLevel))
      # check that the desired contrast is already
      # available in mcols(object), and then we can either
      # take it directly or multiply the log fold
      # changes and Wald stat by -1
      if ( contrastDenomLevel == contrastBaseLevel ) {
        cleanName <- paste(contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
        # the results can be pulled directly from mcols(object)
        name <- if (!noInterceptPullCoef) {
          make.names(paste0(contrastFactor,"_",contrastNumLevel,"_vs_",contrastDenomLevel))
        } else {
          make.names(paste0(contrastFactor,contrastNumLevel))
        }
        if (!name %in% resNames) {
          stop(paste("as",contrastDenomLevel,"is the reference level, was expecting",name,"to be present in 'resultsNames(object)'"))
        }
        log2FoldChange <- getCoef(object, name)
        lfcSE <- getCoefSE(object, name)
        stat <- getStat(object, test, name)
        pvalue <- getPvalue(object, test, name)
        res <- cbind(mcols(object)["baseMean"],log2FoldChange,lfcSE,stat,pvalue)
        names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
        lfcType <- if (attr(object,"betaPrior")) "MAP" else "MLE"
        lfcDesc <- paste0("log2 fold change (",lfcType,"): ",cleanName)
        mcols(res,use.names=TRUE)["log2FoldChange","description"] <- lfcDesc
        resReady <- TRUE
        
      } else if ( contrastNumLevel == contrastBaseLevel ) {
        # fetch the results for denom vs num 
        # and mutiply the log fold change and stat by -1
        cleanName <- paste(contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
        swapName <- if (!noInterceptPullCoef) {
          make.names(paste0(contrastFactor,"_",contrastDenomLevel,"_vs_",contrastNumLevel))
        } else {
          make.names(paste0(contrastFactor,contrastDenomLevel))
        }
        if (!swapName %in% resNames) {
          stop(paste("as",contrastNumLevel,"is the reference level, was expecting",swapName,"to be present in 'resultsNames(object)'"))
        }
        log2FoldChange <- getCoef(object, swapName)
        lfcSE <- getCoefSE(object, swapName)
        stat <- getStat(object, test, swapName)
        pvalue <- getPvalue(object, test, swapName)
        res <- cbind(mcols(object)["baseMean"],log2FoldChange,lfcSE,stat,pvalue)
        names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
        res$log2FoldChange <- -1 * res$log2FoldChange
        if (test == "Wald") res$stat <- -1 * res$stat
        lfcType <- if (attr(object,"betaPrior")) "MAP" else "MLE"
        # rename some of the columns using the flipped contrast
        if (test == "Wald") {
          contrastDescriptions <- paste(c(paste0("log2 fold change (",lfcType,"):"),
                                          "standard error:",
                                          "Wald statistic:","Wald test p-value:"), cleanName)
          mcols(res,use.names=TRUE)[c("log2FoldChange","lfcSE","stat","pvalue"),
                      "description"] <- contrastDescriptions
        } else {
          contrastDescriptions <- paste(c(paste0("log2 fold change (",lfcType,"):"),
                                          "standard error:"), cleanName)
          mcols(res,use.names=TRUE)[c("log2FoldChange","lfcSE"),
                      "description"] <- contrastDescriptions
        }
        resReady <- TRUE
        
      } else {
        # check for the case where neither are present
        # as comparisons against reference level
        if ( ! (contrastNumColumn %in% resNames &
                  contrastDenomColumn %in% resNames) ) {
          stop(paste(contrastNumLevel,"and",contrastDenomLevel,"should be levels of",contrastFactor,"such that",contrastNumColumn,"and",contrastDenomColumn,"are contained in 'resultsNames(object)'"))
        }
      }
      # case 2: expanded model matrices or no intercept and first variable
      # need to then build the appropriate contrast.
      # these coefficient names have the form "factorlevel"
      # output: contrastNumColumn and contrastDenomColumn
    } else {
      # we only need to check validity
      contrastNumColumn <- make.names(paste0(contrastFactor, contrastNumLevel))
      contrastDenomColumn <- make.names(paste0(contrastFactor, contrastDenomLevel))
      if ( ! (contrastNumColumn %in% resNames & contrastDenomColumn %in% resNames) ) {
        stop(paste(paste0(contrastFactor,contrastNumLevel),"and",paste0(contrastFactor,contrastDenomLevel),
                   "are expected to be in resultsNames(object)"))
      }
    }

    # check if both levels have all zero counts
    # (this has to be down here to make use of error checking above)
    contrastAllZero <- contrastAllZeroCharacter(object, contrastFactor,
                                                contrastNumLevel, contrastDenomLevel)
    
  }

  # if the result table not already built in the above code
  if (!resReady) {
    
    # here, a numeric / list / character contrast which will be converted
    # into a numeric contrast and run through getContrast()
    if (is.numeric(contrast)) {
      # make name for numeric contrast
      signMap <- c("","","+")
      contrastSigns <- signMap[sign(contrast)+2]
      contrastName <- paste(paste0(contrastSigns,as.character(contrast)),collapse=",")
    } else if (is.list(contrast)) {
      # interpret list contrast into numeric and make a name for the contrast
      lc1 <- length(contrast[[1]])
      lc2 <- length(contrast[[2]])
      # these just used for naming
      listvalname1 <- round(listValues[1],3)
      listvalname2 <- round(listValues[2],3)
      if (lc1 > 0 & lc2 > 0) {
        listvalname2 <- abs(listvalname2)
        listvalname1 <- if (listvalname1 == 1) "" else paste0(listvalname1," ")
        listvalname2 <- if (listvalname2 == 1) "" else paste0(listvalname2," ")
        contrastName <- paste0(listvalname1,paste(contrast[[1]],collapse="+")," vs ",listvalname2,paste(contrast[[2]],collapse="+"))
      } else if (lc1 > 0 & lc2 == 0) {
        listvalname1 <- if (listvalname1 == 1) "" else paste0(listvalname1," ")
        contrastName <- paste0(listvalname1,paste(contrast[[1]],collapse="+")," effect")
      } else if (lc1 == 0 & lc2 > 0) {
        contrastName <- paste(listvalname2,paste(contrast[[2]],collapse="+"),"effect")
      }
      contrastNumeric <- rep(0,length(resNames))
      contrastNumeric[resNames %in% contrast[[1]]] <- listValues[1]
      contrastNumeric[resNames %in% contrast[[2]]] <- listValues[2]
      contrast <- contrastNumeric
    } else if (is.character(contrast)) {
      # interpret character contrast into numeric and make a name for the contrast
      contrastNumeric <- rep(0,length(resNames))
      contrastNumeric[resNames == contrastNumColumn] <- 1
      contrastNumeric[resNames == contrastDenomColumn] <- -1
      contrast <- contrastNumeric
      contrastName <- paste(contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
    }

    contrastAllZero <- contrastAllZeroNumeric(object, contrast)

    # now get the contrast
    contrastResults <- getContrast(object, contrast, useT=useT, minmu)
    lfcType <- if (attr(object,"betaPrior")) "MAP" else "MLE"
    contrastDescriptions <- paste(c(paste0("log2 fold change (",lfcType,"):"),
                                    "standard error:",
                                    "Wald statistic:",
                                    "Wald test p-value:"),
                                  contrastName)
    mcols(contrastResults) <- DataFrame(type=rep("results",ncol(contrastResults)),
                                        description=contrastDescriptions)
    res <- cbind(mcols(object)["baseMean"],
                 contrastResults)

  }
    
  # if the counts in all samples included in contrast are zero
  # then zero out the LFC, Wald stat and p-value set to 1
  contrastAllZero <- contrastAllZero & !mcols(object)$allZero
  if (sum(contrastAllZero) > 0) {
    res$log2FoldChange[contrastAllZero] <- 0
    res$stat[contrastAllZero] <- 0
    res$pvalue[contrastAllZero] <- 1
  }

  # if test is "LRT", overwrite the statistic and p-value
  # (we only ran contrast for the coefficient)
  if (test == "LRT") {
    stat <- getStat(object, test, name=NULL)
    pvalue <- getPvalue(object, test, name=NULL)
    res <- cbind(res[c("baseMean","log2FoldChange","lfcSE")],stat,pvalue)
    names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
  }
  
  return(res)
}


# convenience function to guess the name of the last coefficient
# in the model matrix, unless specified this will be used for
# plots and accessor functions
lastCoefName <- function(object) {
  resNms <- resultsNames(object)
  resNms[length(resNms)]
}

# functions to get coef, coefSE, pvalues and padj from mcols(object)
getCoef <- function(object,name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  mcols(object)[name]
}
getCoefSE <- function(object,name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  mcols(object)[paste0("SE_",name)]
}
getStat <- function(object,test="Wald",name) {
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  if (test == "Wald") {
    return(mcols(object)[paste0("WaldStatistic_",name)])
  } else if (test == "LRT") {
    return(mcols(object)["LRTStatistic"])
  } else {
    stop("unknown test")
  }
}
getPvalue <- function(object,test="Wald",name) {
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

# convenience function to make more descriptive names
# for factor variables
renameModelMatrixColumns <- function(data, design) {
  data <- as.data.frame(data)
  designVars <- all.vars(design)
  designVarsClass <- sapply(designVars, function(v) class(data[[v]]))
  factorVars <- designVars[designVarsClass == "factor"]
  colNamesFrom <- make.names(do.call(c,lapply(factorVars, function(v) paste0(v,levels(data[[v]])[-1]))))
  colNamesTo <- make.names(do.call(c,lapply(factorVars, function(v) paste0(v,"_",levels(data[[v]])[-1],"_vs_",levels(data[[v]])[1]))))
  data.frame(from=colNamesFrom,to=colNamesTo,stringsAsFactors=FALSE)
}

makeWaldTest <- function(object) {
  betaMatrix <- as.matrix(mcols(object)[,grep("log2 fold change",mcols(mcols(object))$description),drop=FALSE])
  modelMatrixNames <- colnames(betaMatrix)
  betaSE <- as.matrix(mcols(object)[,grep("standard error",mcols(mcols(object))$description),drop=FALSE])
  WaldStatistic <- betaMatrix/betaSE
  colnames(WaldStatistic) <- paste0("WaldStatistic_",modelMatrixNames)
  # don't need to worry about t distribution here bc nbniomLRT() was run
  WaldPvalue <- 2*pnorm(abs(WaldStatistic),lower.tail=FALSE)
  colnames(WaldPvalue) <- paste0("WaldPvalue_",modelMatrixNames)
  modelMatrixNamesSpaces <- gsub("_"," ",modelMatrixNames)
  statInfo <- paste("Wald statistic:",modelMatrixNamesSpaces)
  pvalInfo <- paste("Wald test p-value:",modelMatrixNamesSpaces)
  WaldResults <- DataFrame(c(matrixToList(WaldStatistic), matrixToList(WaldPvalue)))
  mcols(WaldResults) <- DataFrame(type = rep("results",ncol(WaldResults)),
                                  description = c(statInfo, pvalInfo))
  mcols(object) <- cbind(mcols(object),WaldResults)
  return(object)
}

mleContrast <- function(object, contrast) {
  contrastFactor <- contrast[1]
  contrastNumLevel <- contrast[2]
  contrastDenomLevel <- contrast[3]
  contrastRefLevel <- levels(colData(object)[,contrastFactor])[1]
  contrastNumColumn <- make.names(paste0("MLE_",contrastFactor,"_",contrastNumLevel,"_vs_",contrastRefLevel))
  contrastDenomColumn <- make.names(paste0("MLE_",contrastFactor,"_",contrastDenomLevel,"_vs_",contrastRefLevel))
  cleanName <- paste("log2 fold change (MLE):",contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
  if ( contrastDenomLevel == contrastRefLevel ) {
    name <- make.names(paste0("MLE_",contrastFactor,"_",contrastNumLevel,"_vs_",contrastDenomLevel))
    lfcMLE <- mcols(object)[name]
  } else if ( contrastNumLevel == contrastRefLevel ) {
    swapName <- make.names(paste0("MLE_",contrastFactor,"_",contrastDenomLevel,"_vs_",contrastNumLevel))
    lfcMLE <- mcols(object)[swapName]
    lfcMLE[[1]] <- -1 * lfcMLE[[swapName]]
  } else {
    numMLE <- mcols(object)[[contrastNumColumn]]
    denomMLE <- mcols(object)[[contrastDenomColumn]]
    lfcMLE <- mcols(object)[contrastNumColumn]
    lfcMLE[[1]] <- numMLE - denomMLE
  }
  names(lfcMLE) <- "lfcMLE"
  mcols(lfcMLE)$description <- cleanName
  lfcMLE
}

checkContrast <- function(contrast, resNames) {
  if (!(is.numeric(contrast) | is.character(contrast) | is.list(contrast))) {
    stop("'contrast' vector should be either a character vector of length 3,
a list of length 2 containing character vectors,
or a numeric vector, see the argument description in ?results")
  }

  # character
  if (is.character(contrast)) {
    if (length(contrast) != 3) {
      stop("'contrast', as a character vector of length 3, should have the form:
contrast = c('factorName','numeratorLevel','denominatorLevel'),
see the manual page of ?results for more information")
    }
    if (contrast[2] == contrast[3]) {
      stop(paste(contrast[2],"and",contrast[3],"should be different level names"))
    }
  }

  # list
  if (is.list(contrast)) {
    if (length(contrast) == 1) {
      contrast <- list(contrast[[1]], character())
    }
    if (length(contrast) != 2) {
      stop("'contrast', as a list, should have length 2, or, if length 1,
an empty vector will be added for the second element.
see the manual page of ?results for more information")
    }
    if (!(is.character(contrast[[1]]) & is.character(contrast[[2]]))) {
      stop("'contrast', as a list of length 2, should have character vectors as elements,
see the manual page of ?results for more information")
    }
    if (!all(c(contrast[[1]],contrast[[2]]) %in% resNames)) {
      stop("all elements of the contrast as a list of length 2 should be elements of 'resultsNames(object)'")
    }
    if (length(intersect(contrast[[1]], contrast[[2]])) > 0) {
      stop("elements in the contrast list should only appear in the numerator (first element of contrast list)
or the denominator (second element of contrast list), but not both")
    }
    if (length(c(contrast[[1]],contrast[[2]])) == 0) {
      stop("one of the two elements in the list should be a character vector of non-zero length")
    }    
  }

  # numeric
  if (is.numeric(contrast)) {
    if (length(contrast) != length(resNames) )
      stop("numeric contrast vector should have one element for every element of 'resultsNames(object)'")
    if (all(contrast==0)) {
      stop("numeric contrast vector cannot have all elements equal to 0")
    }
  }

  return(contrast)
}

# function to determine output of results() and lfcShrink()
resultsFormatSwitch <- function(object, res, format, saveCols) {
  
  if (format == "DataFrame") {
    out <- res
  } else if (format == "GRangesList") {
    if (is(rowRanges(object), "GRanges")) warning("rowRanges is GRanges")
    out <- rowRanges(object)
    mcols(out) <- res
  } else if (format == "GRanges") {
    if (is(rowRanges(object), "GRangesList")) {
      message("rowRanges is GRangesList, performing unlist(range(x)) on the rowRanges")
      out <- unlist(range(rowRanges(object)))
      mcols(out) <- res
    } else {
      out <- rowRanges(object)
      mcols(out) <- res
    }
  }

  # code to save certain columns of mcols(rowRanges) / rowData
  if (!is.null(saveCols)) {
    mcols2Save <- mcols(object)[,saveCols,drop=FALSE]
    if (format == "DataFrame") {
      out <- cbind(out, mcols2Save)
    } else {
      mcols(out) <- cbind(mcols(out), mcols2Save)
    }
  }

  out
}


contrastAllZeroCharacter <- function(object, contrastFactor,
                                     contrastNumLevel, contrastDenomLevel) {
  cts <- counts(object)
  f <- colData(object)[[contrastFactor]]
  cts.sub <- cts[ , f %in% c(contrastNumLevel, contrastDenomLevel), drop=FALSE ]
  rowSums( cts.sub == 0 ) == ncol(cts.sub)
}

contrastAllZeroNumeric <- function(object, contrast) {
  if (is.null(attr(object,"modelMatrix"))) {
    stop("was expecting a model matrix stored as an attribute of the DESeqDataSet")
  }
  modelMatrix <- attr(object, "modelMatrix")
  # note: this extra leg-work to zero out LFC, lfcSE, and set p-value to 1
  # for contrasts comparing groups where both groups have all zeros.

  # it is only implemented for the case in which we can identify
  # the relevant samples by multiplying the model matrix
  # with a vector where the non-zero elements of the numeric contrast
  # are replaced with 1.

  # so, this code will not zero out in the case of standard model matrices
  # where the user supplies a numeric contrast that pulls out a single column
  # of the model matrix, for example.
  
  if (all(contrast >= 0) | all(contrast <= 0)) {
    return( rep(FALSE, nrow(object)) )
  }
  
  contrastBinary <- ifelse(contrast == 0, 0, 1)
  whichSamples <- ifelse(modelMatrix %*% contrastBinary == 0, 0, 1)
  zeroTest <- counts(object) %*% whichSamples
  zeroTest == 0
}
