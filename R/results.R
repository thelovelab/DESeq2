#' Extract results from a DESeq analysis
#'
#' \code{results} extracts results from a DESeq analysis
#' giving base means across samples,
#' log2 fold changes, standard errors, test statistics,
#' p-values and adjusted p-values.
#' \code{resultsNames} finds available names for results;
#' \code{removeResults} returns a \code{DESeqDataSet} object
#' with results columns removed.
#'
#' Multiple results can be returned for an analysis with multifactor design,
#' so \code{results} takes an argument \code{name} as well as the
#' argument \code{contrast} explained below.
#'
#' The available names can be checked using \code{resultsNames};
#' these are combined variable names and factor levels, potentially 
#' with minor changes made by the \code{make.names} function on column names
#' (e.g. dashes into periods).
#'
#' By default, results for the last variable will be returned.
#' Information on the variable represented and the test
#' used for p-values (Wald test or likelihood ratio test)
#' is stored in the metadata columns, accessible by calling \code{mcol}
#' on the \code{DataFrame} returned by \code{results}.
#'
#' By default, independent filtering is performed to select a set of genes
#' which will result in the most genes with adjusted p-values less than a
#' threshold, alpha. The adjusted p-values for the genes which do not pass 
#' the filter threshold are set to \code{NA}. By default, the mean of normalized counts
#' is used to perform this filtering, though other statistics can be provided.
#' Several arguments from the \code{filtered_p} function of genefilter
#' are provided to control or turn off the independent filtering behavior.
#'
#' A contrast can be performed by specifying the variable and factor
#' levels which should be compared, or by specifying the numeric contrast
#' vector. The test statistic is then:
#'
#' \deqn{ c^t \beta / \sqrt{c^t \Sigma c } }{ c' beta / sqrt( c' Sigma c ) }
#'
#' For analyses using the likelihood ratio test (using \code{\link{nbinomLRT}}),
#' the p-values are determined solely by the difference in deviance between
#' the full and reduced model formula.
#' In this case, the \code{name} argument only specifies which
#' coefficient should be used for reporting the log2 fold changes.
#'
#' Cook's distance for each sample are accessible as a matrix "cooks"
#' stored in the assays() list. This measure is useful for identifying rows where the
#' observed counts might not fit to a negative binomial distribution.
#' 
#' Results can be removed from an object by calling \code{removeResults}
#'
#' @references Richard Bourgon, Robert Gentleman, Wolfgang Huber: Independent
#' filtering increases detection power for high-throughput experiments.
#' PNAS (2010), \url{http://dx.doi.org/10.1073/pnas.0914005107}
#' 
#' @param object a DESeqDataSet, on which one
#' of the following functions has already been called:
#' \code{\link{DESeq}}, \code{\link{nbinomWaldTest}}, or \code{\link{nbinomLRT}}
#' @param name the name of the coefficient for which to report log2 fold changes
#' -- and for the Wald test, p-values and adjusted p-values
#' @param contrast either a character vector with exactly three elements
#' (name of the factor, name of the numerator level, name of the
#' denominator level), or a numeric contrast vector with one element
#' for each element in \code{resultsNames(object)}, i.e. columns of the model matrix.
#' The DESeqDataSet must be one produced using the Wald test steps
#' in order to use contrasts. If used, the \code{name} argument is ignored.
#' @param lfcThreshold a non-negative value, which specifies the test which should
#' be applied to the log2 fold changes. The standard is a test that the log2 fold
#' changes are not equal to zero. However, log2 fold changes greater or less than
#' \code{lfcThreshold} can also be tested. Specify the alternative hypothesis
#' using the \code{altHypothesis} argument. If \code{lfcThreshold} is specified,
#' the results are Wald tests, and LRT p-values will be overwritten.
#' @param altHypothesis character which specifies the alternative hypothesis,
#' i.e. those values of log2 fold change which the user is interested in
#' finding. The complement of this set of values is the null hypothesis which
#' will be tested. If the log2 fold change specified by \code{name}
#' or by \code{contrast} is written as \eqn{ \beta }{ beta }, then the possible values for
#' \code{altHypothesis} represent the following alternate hypotheses:
#' \itemize{
#' \item greaterAbs - \eqn{|\beta| > \textrm{lfcThreshold} }{ |beta| > lfcThreshold },
#' and p-values are two-tailed
#' \item lessAbs - \eqn{ |\beta| < \textrm{lfcThreshold} }{ |beta| < lfcThreshold },
#' NOTE: this requires that \code{betaPrior=FALSE} has been specified in the 
#' previous \code{\link{DESeq}} call. 
#' p-values are the maximum of the upper and lower tests.
#' \item greater - \eqn{ \beta > \textrm{lfcThreshold} }{ beta > lfcThreshold }
#' \item less - \eqn{ \beta < -\textrm{lfcThreshold} }{ beta < -lfcThreshold }
#' }
#' @param cooksCutoff theshold on Cook's distance, such that if one or more
#' samples for a row have a distance higher, the p-value for the row is
#' set to NA.
#' The default cutoff is the .99 quantile of the F(p, m-p) distribution,
#' where p is the number of coefficients being fitted and m is the number of samples.
#' Set to Inf or FALSE to disable the resetting of p-values to NA.
#' Note: this test excludes the Cook's distance of samples whose removal
#' would result in rank deficient design matrix and samples belonging to experimental
#' groups with only 2 samples.
#' @param independentFiltering logical, whether independent filtering should be
#' applied automatically
#' @param alpha the significance cutoff used for optimizing the independent
#' filtering
#' @param filter the vector of filter statistics over which the independent
#' filtering will be optimized. By default the mean of normalized counts is used.
#' @param theta the quantiles at which to assess the number of rejections
#' from independent filtering
#' @param pAdjustMethod the method to use for adjusting p-values, see \code{?p.adjust}
#'
#' @return For \code{results}: a DataFrame of results columns:
#' \code{log2FoldChange}, \code{lfcSE}, \code{stat},
#' \code{pvalue} and \code{padj}.
#' also includes metadata columns of variable information
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
results <- function(object, name, contrast,
                    lfcThreshold=0,
                    altHypothesis=c("greaterAbs","lessAbs","greater","less"),
                    cooksCutoff,
                    independentFiltering=TRUE,
                    alpha=0.1, filter, theta,
                    pAdjustMethod="BH") {
  if (!"results" %in% mcols(mcols(object))$type) {
    stop("cannot find results columns in object, first call 'DESeq','nbinomWaldTest', or 'nbinomLRT'")
  }
  isExpanded <- attr(object, "modelMatrixType") == "expanded"
  termsOrder <- attr(terms.formula(design(object)),"order")
  if (isExpanded & missing(contrast) & all(termsOrder < 2)) {
    if (missing(name)) {
      designVars <- all.vars(formula(design(object)))
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      nlvls <- nlevels(lastVar)
      contrast <- c(lastVarName, levels(lastVar)[nlvls], levels(lastVar)[1])
    } else {     
      message("\n
note: an expanded model matrix was used in fitting, either through
use of the modelMatrixType argument or by default, because 1 or
more factors in the design formula contained 3 or more levels.

recommendation: the 'contrast' argument should be used to extract
log2 fold changes of levels against each other, otherwise the log2
fold changes are compared to the intercept.\n")
    }
  }
  if (missing(name)) {
    name <- lastCoefName(object)
  }
  altHypothesis <- match.arg(altHypothesis, choices=c("greaterAbs","lessAbs","greater","less"))
  stopifnot(lfcThreshold >= 0)
  stopifnot(length(lfcThreshold)==1)
  stopifnot(length(altHypothesis)==1)
  stopifnot(length(alpha)==1)
  stopifnot(length(pAdjustMethod)==1)
  if (length(name) != 1 | !is.character(name)) {
    stop("the argument 'name' should be a character vector of length 1")
  }
  if (lfcThreshold == 0 & altHypothesis == "lessAbs") {
    stop("when testing altHypothesis='lessAbs', set the argument lfcThreshold to a positive value")
  }
  
  # determine test type from the names of mcols(object)
  if (paste0("WaldPvalue_",name) %in% names(mcols(object))) {
    test <- "Wald"
  } else if ("LRTPvalue" %in% names(mcols(object))) {
    test <- "LRT"
  } else {
    stop("cannot find appropriate results, for available names call 'resultsNames(object)'")
  }
  
  # if performing a contrast call the function cleanContrast()
  if (!missing(contrast)) {
    # must have performed the Wald test steps
    if (test != "Wald") {
      stop("using contrasts requires that the Wald test was performed")
    }
    # pass down whether the model matrix type was "expanded"
    res <- cleanContrast(object, contrast, expanded=isExpanded)
  } else {
    # if not performing a contrast
    # pull relevant columns from mcols(object)
    log2FoldChange <- getCoef(object, name)
    lfcSE <- getCoefSE(object, name)
    stat <- getStat(object, test, name)
    pvalue <- getPvalue(object, test, name)
    res <- cbind(mcols(object)["baseMean"],
                 log2FoldChange,lfcSE,stat,
                 pvalue)
    names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
  }
  rownames(res) <- rownames(object)
  
  # only if we need to generate new p-values
  if ( !(lfcThreshold == 0 & altHypothesis == "greaterAbs") ) {
    if (test == "LRT") {
      warning("tests of log fold change above or below a theshold are Wald tests. Likelihood ratio test p-values are overwritten")
    }
    if (altHypothesis == "greaterAbs") {
      newPvalue <- 2 * pnorm(abs(res$log2FoldChange),
                             mean = lfcThreshold,
                             sd = res$lfcSE,
                             lower.tail = FALSE)
    } else if (altHypothesis == "lessAbs") {
      # check requirement if betaPrior was set to FALSE
      if (attr(object,"betaPrior")) {
        stop("testing altHypothesis='lessAbs' requires setting the DESeq() argument betaPrior=FALSE")
      }
      pvalueAbove <- pnorm(res$log2FoldChange,
                           mean = lfcThreshold,
                           sd = res$lfcSE,
                           lower.tail = TRUE)
      pvalueBelow <- pnorm(res$log2FoldChange,
                           mean = -lfcThreshold,
                           sd = res$lfcSE,
                           lower.tail = FALSE)
      newPvalue <- pmax(pvalueAbove, pvalueBelow)
    } else if (altHypothesis == "greater") {
      newPvalue <- pnorm(res$log2FoldChange,
                         mean = lfcThreshold,
                         sd = res$lfcSE,
                         lower.tail = FALSE)
    } else if (altHypothesis == "less") {
      newPvalue <- pnorm(res$log2FoldChange,
                         mean = -lfcThreshold,
                         sd = res$lfcSE,
                         lower.tail = TRUE)
    }
    res$pvalue <- newPvalue
  }
  
  # calculate Cook's cutoff
  m <- nrow(attr(object,"modelMatrix"))
  p <- ncol(attr(object,"modelMatrix"))
  
  # only if more samples than parameters:
  if (m > p) {
    defaultCutoff <- qf(.99, p, m - p)
    if (missing(cooksCutoff)) {
      cooksCutoff <- defaultCutoff
    }
    stopifnot(length(cooksCutoff)==1)
    if (is.logical(cooksCutoff) & cooksCutoff) {
      cooksCutoff <- defaultCutoff
    }
  } else {
    cooksCutoff <- FALSE
  }
  
  # apply cutoff based on maximum Cook's distance
  performCooksCutoff <- (is.numeric(cooksCutoff) | cooksCutoff) 
  if ((m > p) & performCooksCutoff) {
    cooksOutlier <- mcols(object)$maxCooks > cooksCutoff
    res$pvalue[cooksOutlier] <- NA
  }
  
  # perform independent filtering
  if (independentFiltering) {
    if (missing(filter)) {
      filter <- res$baseMean
    }
    if (missing(theta)) {
      lowerQuantile <- mean(filter == 0)
      if (lowerQuantile < .95) upperQuantile <- .95 else upperQuantile <- 1
      theta <- seq(lowerQuantile, upperQuantile, length=20)
    }
    stopifnot(length(theta) > 1)
    stopifnot(length(filter) == nrow(object))
    filtPadj <- filtered_p(filter=filter, test=res$pvalue,
                           theta=theta, method=pAdjustMethod) 
    numRej  <- colSums(filtPadj < alpha, na.rm = TRUE)
    j <- which.max(numRej)
    res$padj <- filtPadj[, j, drop=TRUE]
    cutoffs <- quantile(filter, theta) 
    attr(res, "filterThreshold") <- cutoffs[j]
    attr(res, "filterNumRej") <- data.frame(theta=theta, numRej=numRej)
  } else {
    # regular p-value adjustment
    # which does not include those rows which were removed
    # by maximum Cook's distance
    res$padj <- p.adjust(res$pvalue,method=pAdjustMethod)
  }
  
  mcols(res)$type[names(res)=="padj"] <- "results"
  mcols(res)$description[names(res)=="padj"] <- paste(pAdjustMethod,"adjusted p-values")
  
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
    mcols(object) <- mcols(object)[,-which(resCols)]
  }
  return(object)
}

###########################################################
# unexported functons 
###########################################################

# unexported function
# results() just calls filtered_p directly
pAdjustWithIndependentFiltering <- function(p, filterstat, alpha=0.1, method = "BH", plot=FALSE ){
  theta <- seq(0, 0.95, by=0.05)
  cutoffs <- quantile( filterstat, theta ) 
  padj <- filtered_p( filterstat, p, cutoffs, method ) 
  rej  <- colSums(padj < alpha, na.rm = TRUE )
  if(plot) plot(theta, rej)
  j <- which.max(rej)
  rv <- padj[, j, drop=TRUE]
  attr(rv, "filterthreshold") = cutoffs[j]
  rv
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
getContrast <- function(object, contrast, useT=FALSE, df) {
  if (missing(contrast)) {
    stop("must provide a contrast")
  }

  if (is.null(attr(object,"modelMatrix"))) {
    stop("was expecting a model matrix stored as an attribute of the DESeqDataSet")
  }
  modelMatrix <- attr(object, "modelMatrix")
  
  # only continue on the rows with non-zero row mean
  objectNZ <- object[!mcols(object)$allZero,]
  if (!is.null(normalizationFactors(objectNZ))) {
    normalizationFactors <- normalizationFactors(objectNZ)
  } else { 
    normalizationFactors <- matrix(rep(sizeFactors(objectNZ),each=nrow(objectNZ)),
                                   ncol=ncol(objectNZ))
  }
  alpha_hat <- dispersions(objectNZ)
  coefColumns <- names(mcols(objectNZ))[grep("log2 fold change",mcols(mcols(object))$description)]
  # convert betas to log scale
  beta_mat <- log(2) * as.matrix(mcols(objectNZ)[,coefColumns,drop=FALSE])
  # convert beta prior variance to log scale
  lambda = 1/(log(2)^2 * attr(object,"betaPriorVar")) 
  betaRes <- fitBeta(ySEXP = counts(objectNZ), xSEXP = modelMatrix,
                     nfSEXP = normalizationFactors,
                     alpha_hatSEXP = alpha_hat,
                     contrastSEXP = contrast,
                     beta_matSEXP = beta_mat,
                     lambdaSEXP = lambda,
                     tolSEXP = 1e-8, maxitSEXP = 0,
                     useQRSEXP=FALSE)
  # convert back to log2 scale
  contrastEstimate <- log2(exp(1)) * betaRes$contrast_num
  contrastSE <- log2(exp(1)) * betaRes$contrast_denom
  contrastStatistic <- contrastEstimate / contrastSE
  if (useT) {
    stopifnot(length(df)==1)
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

cleanContrast <- function(object, contrast, expanded=FALSE) {
  # get the names of columns in the beta matrix
  resNames <- resultsNames(object)
  
  if (!is.numeric(contrast) & !is.character(contrast)) {
    stop("contrast vector should be either a numeric vector or character vector,
see the argument description in ?results")
  }
  
  # check contrast validity
  if (is.numeric(contrast)) {
    if (length(contrast) != length(resNames) )
      stop("numeric contrast vector should have one element for every element of 'resultsNames(object)'")
    if (all(contrast==0)) {
      stop("numeric contrast vector cannot have all elements equal to 0")
    }
  } else if (is.character(contrast)) {
    # check if the appropriate columns are in the resultsNames
    if (contrast[2] == contrast[3]) {
      stop(paste(contrast[2],"and",contrast[3],"should be different level names"))
    }
    contrastFactor <- contrast[1]
    if (!contrastFactor %in% names(colData(object))) {
      stop(paste(contrastFactor,"should be the name of a factor in the colData of the DESeqDataSet"))
    }
    contrastNumLevel <- contrast[2]
    contrastDenomLevel <- contrast[3]

    # more checks, and potentially return results already
    # if numerator or denominator is the base level
    # first for standard model matrices
    if (!expanded) {

      # then we have a base level for the factor
      contrastBaseLevel <- levels(colData(object)[,contrastFactor])[1]
      
      # use make.names() so the column names are
      # the same as created by DataFrame in mcols(object).
      contrastNumColumn <- make.names(paste0(contrastFactor,"_",contrastNumLevel,"_vs_",contrastBaseLevel))
      contrastDenomColumn <- make.names(paste0(contrastFactor,"_",contrastDenomLevel,"_vs_",contrastBaseLevel))
      resNames <- resultsNames(object)
      
      # check in case the desired contrast is already
      # available in mcols(object), and then we can either
      # take it directly or multiply the log fold
      # changes and stat by -1
      if ( contrastDenomLevel == contrastBaseLevel ) {
        # the results can be pulled directly from mcols(object)
        name <- make.names(paste0(contrastFactor,"_",contrastNumLevel,"_vs_",contrastDenomLevel))
        if (!name %in% resNames) {
          stop(paste("as",contrastDenomLevel,"is the base level, was expecting",name,"to be present in 'resultsNames(object)'"))
        }
        test <- "Wald"
        log2FoldChange <- getCoef(object, name)
        lfcSE <- getCoefSE(object, name)
        stat <- getStat(object, test, name)
        pvalue <- getPvalue(object, test, name)
        res <- cbind(mcols(object)["baseMean"],
                     log2FoldChange,lfcSE,stat,
                     pvalue)
        names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
        return(res)
      } else if ( contrastNumLevel == contrastBaseLevel ) {
        # fetch the results for denom vs num 
        # and mutiply the log fold change and stat by -1
        cleanName <- paste(contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
        swapName <- make.names(paste0(contrastFactor,"_",contrastDenomLevel,"_vs_",contrastNumLevel))
        if (!swapName %in% resNames) {
          stop(paste("as",contrastNumLevel,"is the base level, was expecting",swapName,"to be present in 'resultsNames(object)'"))
        }
        test <- "Wald"
        log2FoldChange <- getCoef(object, swapName)
        lfcSE <- getCoefSE(object, swapName)
        stat <- getStat(object, test, swapName)
        pvalue <- getPvalue(object, test, swapName)
        res <- cbind(mcols(object)["baseMean"],
                     log2FoldChange,lfcSE,stat,
                     pvalue)
        names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
        res$log2FoldChange <- -1 * res$log2FoldChange
        res$stat <- -1 * res$stat
        # also need to swap the name in the mcols
        contrastDescriptions <- paste(c("log2 fold change (MAP):",
                                        "standard error:",
                                        "Wald statistic:",
                                        "Wald test p-value:"),
                                      cleanName)
        mcols(res)$description[mcols(res)$type == "results"] <- contrastDescriptions
        return(res)
      }
    
      # check for the case where neither are present
      # as comparisons against base level
      if ( ! (contrastNumColumn %in% resNames &
              contrastDenomColumn %in% resNames) ) {
        stop(paste(contrastNumLevel,"and",contrastDenomLevel,"should be levels of",contrastFactor,
                   "such that",contrastNumColumn,"and",contrastDenomColumn,
                   "are contained in 'resultsNames(object)'"))
      }
      # then we proceed below with translation to a numeric contrast
      
    } else {
      
      # else in the expanded case, we first check validity
      contrastNumColumn <- make.names(paste0(contrastFactor, contrastNumLevel))
      contrastDenomColumn <- make.names(paste0(contrastFactor, contrastDenomLevel))
      if ( ! (contrastNumColumn %in% resNames & contrastDenomColumn %in% resNames) ) {
        stop(paste("both",contrastNumLevel,"and",contrastDenomLevel,"are expected to be in
resultsNames(object), prefixed by",contrastFactor))
      }
    }
  }

  # make name for numeric contrast
  if (is.numeric(contrast)) {
    signMap <- c("","","+")
    contrastSigns <- signMap[sign(contrast)+2]
    contrastName <- paste(paste0(contrastSigns,as.character(contrast)),collapse=",")
  } else if (is.character(contrast)) {
    # interpret character contrast into numeric
    # and make a name for the contrast
    contrastNumeric <- rep(0,length(resNames))
    contrastNumeric[resNames == contrastNumColumn] <- 1
    contrastNumeric[resNames == contrastDenomColumn] <- -1
    contrast <- contrastNumeric
    contrastName <- paste(contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
  }
  
  contrastResults <- getContrast(object, contrast, useT=FALSE, df)
  contrastDescriptions <- paste(c("log2 fold change (MAP):",
                                  "standard error:",
                                  "Wald statistic:",
                                  "Wald test p-value:"),
                                contrastName)
  mcols(contrastResults) <- DataFrame(type=rep("results",ncol(contrastResults)),
                                      description=contrastDescriptions)
  res <- cbind(mcols(object)["baseMean"],
               contrastResults)
  res
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
renameModelMatrixColumns <- function(modelMatrixNames, data, design) {
  designVars <- all.vars(formula(design))
  designVarsClass <- sapply(designVars, function(v) class(data[[v]]))
  factorVars <- designVars[designVarsClass == "factor"]
  colNamesFrom <- do.call(c,lapply(factorVars, function(v) paste0(v,levels(data[[v]])[-1])))
  colNamesTo <- do.call(c,lapply(factorVars, function(v) paste0(v,"_",levels(data[[v]])[-1],"_vs_",levels(data[[v]])[1])))
  data.frame(from=colNamesFrom,to=colNamesTo,stringsAsFactors=FALSE)
}

getNonInteractionColumnIndices <- function(object, modelMatrix) {
  interactions <- which(attr(terms.formula(design(object)),"order") > 1)
  which(attr(modelMatrix,"assign") != interactions)
}
