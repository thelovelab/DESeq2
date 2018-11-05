#' Shrink log2 fold changes
#'
#' Adds shrunken log2 fold changes (LFC) and SE to a
#' results table from \code{DESeq} run without LFC shrinkage.
#' For consistency with \code{results}, the column name \code{lfcSE}
#' is used here although what is returned is a posterior SD.
#' Three shrinkage esimators for LFC are available via \code{type}.
#'
#' As of DESeq2 version 1.18, \code{type="apeglm"} and \code{type="ashr"}
#' are new features, and still under development.
#' Specifying \code{apeglm} passes along DESeq2 MLE log2
#' fold changes and standard errors to the \code{apeglm} function
#' in the apeglm package, and re-estimates posterior LFCs for
#' the coefficient specified by \code{coef}.
#' Specifying \code{ashr} passes along DESeq2 MLE log2
#' fold changes and standard errors to the \code{ash} function
#' in the ashr package, 
#' with arguments \code{mixcompdist="normal"} and \code{method="shrink"}.
#' See vignette for a comparison of shrinkage estimators on an example dataset.
#' For all shrinkage methods, details on the prior is included in
#' \code{priorInfo(res)}, including the \code{fitted_g} mixture for ashr.
#' The integration of shrinkage methods from
#' external packages will likely evolve over time. We will likely incorporate an
#' \code{lfcThreshold} argument which can be passed to apeglm
#' to specify regions of the posterior at an arbitrary threshold.
#'
#' For \code{normal}, and design as a formula, shrinkage cannot be applied
#' to coefficients in a model with interaction terms. For \code{normal}
#' and user-supplied model matrices, shrinkage is only supported via \code{coef}.
#' For \code{normal} with numeric- or list-style contrasts, it is possible to use \code{lfcShrink},
#' but likely easier to use \code{DESeq} with \code{betaPrior=TRUE} followed by \code{results},
#' because the numeric or list should reference the coefficients from the expanded model matrix.
#' These coefficients will be printed to console if 'contrast' is used with \code{normal}.
#' 
#' @param dds a DESeqDataSet object, after running \code{\link{DESeq}}
#' @param coef the name or number of the coefficient (LFC) to shrink,
#' consult \code{resultsNames(dds)} after running \code{DESeq(dds)}.
#' note: only \code{coef} or \code{contrast} can be specified, not both.
#' \code{apeglm} requires use of \code{coef}.
#' @param contrast see argument description in \code{\link{results}}.
#' only \code{coef} or \code{contrast} can be specified, not both.
#' @param res a DESeqResults object. Results table produced by the
#' default pipeline, i.e. \code{DESeq} followed by \code{results}.
#' If not provided, it will be generated internally using \code{coef} or \code{contrast}.
#' For \code{ashr}, if \code{res} is provided, then \code{coef} and \code{contrast} are ignored.
#' @param type \code{"normal"} is the original DESeq2 shrinkage estimator;
#' \code{"apeglm"} is the adaptive t prior shrinkage estimator from the 'apeglm' package;
#' \code{"ashr"} is the adaptive shrinkage estimator from the 'ashr' package,
#' using a fitted mixture of normals prior
#' - see the Stephens (2016) reference below for citation
#' @param lfcThreshold a non-negative value which specifies a log2 fold change
#' threshold (as in \code{\link{results}}). This can be used in conjunction with
#' \code{normal} and \code{apeglm}, where it will produce new p-values or
#' s-values testing whether the LFC is greater in absolute value than the threshold.
#' The s-values returned in combination with \code{apeglm} provide the probability
#' of FSOS events, "false sign or small", among the tests with equal or smaller s-value
#' than a given gene's s-value, where "small" is specified by \code{lfcThreshold}.
#' @param svalue logical, should p-values and adjusted p-values be replaced
#' with s-values when using \code{apeglm} or \code{ashr}. s-values provide the probability
#' of false signs among the tests with equal or smaller s-value than a given given's s-value.
#' See Stephens (2016) reference on s-values.
#' @param returnList logical, should \code{lfcShrink} return a list, where
#' the first element is the results table, and the second element is the
#' output of \code{apeglm} or \code{ashr}
#' @param format same as defined in \code{\link{results}},
#' either \code{"DataFrame"}, \code{"GRanges"}, or \code{"GRangesList"}
#' @param apeAdapt logical, should \code{apeglm} use the MLE estimates of
#' LFC to adapt the prior, or use default or specified \code{prior.control}
#' @param apeMethod what \code{method} to run \code{apeglm}, which can
#' differ in terms of speed
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see same argument of \code{\link{DESeq}}
#' parallelization only used with \code{normal} or \code{apeglm}
#' @param BPPARAM see same argument of \code{\link{DESeq}}
#' @param quiet whether to print messages 
#' @param ... arguments passed to \code{apeglm} and \code{ashr}
#'
#' @references
#'
#' Publications for the following shrinkage estimators:
#' 
#' \code{type="normal"}:
#'
#' Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. \url{https://doi.org/10.1186/s13059-014-0550-8}
#'
#' \code{type="apeglm"}:
#'
#' Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. Bioinformatics. \url{https://doi.org/10.1093/bioinformatics/bty895}
#' 
#' \code{type="ashr"}:
#'
#' Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. \url{https://doi.org/10.1093/biostatistics/kxw041}
#'
#' Related work, the \code{bayesglm} function in the \code{arm} package:
#'
#' Gelman, A., Jakulin, A., Pittau, M.G. and Su, Y.-S. (2009). A Weakly Informative Default Prior Distribution For Logistic And Other Regression Models. The Annals of Applied Statistics 2:4. \url{http://www.stat.columbia.edu/~gelman/research/published/ priors11.pdf}
#' 
#' @return a DESeqResults object with the \code{log2FoldChange} and \code{lfcSE}
#' columns replaced with shrunken LFC and SE.
#' For consistency with \code{results} (and similar to the output of \code{bayesglm})
#' the column name \code{lfcSE} is used here, although what is returned is a posterior SD.
#' For \code{normal} and for \code{apeglm} the estimate is the posterior mode,
#' for \code{ashr} it is the posterior mean.
#' \code{priorInfo(res)} contains information about the shrinkage procedure,
#' relevant to the various methods specified by \code{type}.
#'
#' @export
#' 
#' @examples
#'
#' set.seed(1)
#' dds <- makeExampleDESeqDataSet(n=500,betaSD=1)
#' dds <- DESeq(dds)
#' res <- results(dds)
#'
#' # these are the coefficients from the model
#' # we can specify them using 'coef' by name or number below
#' resultsNames(dds)
#' 
#' res.shr <- lfcShrink(dds=dds, coef=2)
#' res.shr <- lfcShrink(dds=dds, contrast=c("condition","B","A"))
#' res.ape <- lfcShrink(dds=dds, coef=2, type="apeglm")
#' res.ash <- lfcShrink(dds=dds, coef=2, type="ashr")
#' 
lfcShrink <- function(dds, coef, contrast, res,
                      type=c("normal","apeglm","ashr"),
                      lfcThreshold=0,
                      svalue=FALSE,
                      returnList=FALSE,
                      format=c("DataFrame","GRanges","GRangesList"),
                      apeAdapt=TRUE, apeMethod="nbinomCR",
                      parallel=FALSE, BPPARAM=bpparam(),
                      quiet=FALSE, ...) {  

  stopifnot(is(dds, "DESeqDataSet"))
  if (!missing(res)) {
    if (!is(res, "DESeqResults")) stop("res should be a DESeqResults object, for GRanges output use 'format'")
  }
  type <- match.arg(type, choices=c("normal","apeglm","ashr"))
  format <- match.arg(format, choices=c("DataFrame", "GRanges","GRangesList"))
  if (length(resultsNames(dds)) == 0) {
    stop("first run DESeq() before running lfcShrink()")
  }
  if (attr(dds,"betaPrior")) {
    stop("lfcShrink() should be used downstream of DESeq() with betaPrior=FALSE (the default)")
  }
  stopifnot(length(lfcThreshold) == 1 && lfcThreshold >= 0)
  if (!missing(coef)) {
    if (is.numeric(coef)) {
      stopifnot(coef <= length(resultsNames(dds)))
      coefAlpha <- resultsNames(dds)[coef]
      coefNum <- coef
    } else if (is.character(coef)) {
      stopifnot(coef %in% resultsNames(dds))
      coefNum <- which(resultsNames(dds) == coef)
      coefAlpha <- coef
    }
  }
  if (missing(res)) {
    if (!missing(coef)) {
      res <- results(dds, name=coefAlpha)
    } else if (!missing(contrast)) {
      res <- results(dds, contrast=contrast)
    } else {
      stop("one of coef or contrast required if 'res' is missing")
    }
  }
  if (type %in% c("normal","apeglm")) {
    if (is.null(dispersions(dds))) {
      stop("type='normal' and 'apeglm' require dispersion estimates, first call estimateDispersions()")
    }
    stopifnot(all(rownames(dds) == rownames(res)))
    if (parallel) {
      nworkers <- BPPARAM$workers
      parallelIdx <- factor(sort(rep(seq_len(nworkers),length.out=nrow(dds))))
    }
  }
  
  if (type == "normal") {

    ############
    ## normal ##
    ############

    if (!quiet) message("using 'normal' for LFC shrinkage, the Normal prior from Love et al (2014).
additional priors are available via the 'type' argument, see ?lfcShrink for details")
    
    if (is(design(dds), "formula")) {
      if (attr(dds, "modelMatrixType") == "user-supplied") {
        # if 'full' was used, the model matrix should be stored here
        # TODO... better one day to harmonize these two locations:
        # 1) provided by 'full' and stashed in attr(dds, "modelMatrix")
        # 2) design(dds)
        if (!missing(contrast)) {
          stop("user-supplied design matrix supports shrinkage only with 'coef'")
        }
        modelMatrix <- attr(dds, "modelMatrix")
      } else {
        termsOrder <- attr(terms.formula(design(dds)),"order")
        interactionPresent <- any(termsOrder > 1)
        if (interactionPresent) {
          stop("LFC shrinkage type='normal' not implemented for designs with interactions")
        }
        modelMatrix <- NULL
      }
    } else if (is(design(dds), "matrix")) {
      if (!missing(contrast)) {
        stop("user-supplied design matrix supports shrinkage only with 'coef'")
      }
      modelMatrix <- design(dds)
    }
    
    stopifnot(missing(coef) | missing(contrast))
    # find and rename the MLE columns for estimateBetaPriorVar
    betaCols <- grep("log2 fold change \\(MLE\\)", mcols(mcols(dds))$description)
    stopifnot(length(betaCols) > 0)
    if (!any(grepl("MLE_",names(mcols(dds))[betaCols]))) {
      names(mcols(dds))[betaCols] <- paste0("MLE_", names(mcols(dds))[betaCols])
    }
    if (missing(contrast)) {
      modelMatrixType <- "standard"
    } else {
      # contrast, and so using expanded model matrix: run some checks
      modelMatrixType <- "expanded"
      expMM <- makeExpandedModelMatrix(dds)
      resNames <- colnames(expMM)
      # quick and dirty checks so as to avoid running DESeq() before hitting error
      if (is(contrast, "character")) {
        stopifnot(length(contrast) == 3)
        stopifnot(contrast[1] %in% names(colData(dds)))
        stopifnot(is(colData(dds)[[contrast[1]]], "factor"))
        stopifnot(all(contrast[2:3] %in% levels(colData(dds)[[contrast[1]]])))
      } else {
        message("resultsNames from the expanded model matrix to be referenced by 'contrast':")
        message(paste0("'",paste(resNames, collapse="', '"),"'"))
      }
      contrast <- checkContrast(contrast, resNames)
    }
    attr(dds,"modelMatrixType") <- modelMatrixType
    betaPriorVar <- estimateBetaPriorVar(dds, modelMatrix=modelMatrix)
    stopifnot(length(betaPriorVar) > 0)
    # parallel fork
    if (!parallel) {
      dds.shr <- nbinomWaldTest(dds,
                                betaPrior=TRUE,
                                betaPriorVar=betaPriorVar,
                                modelMatrix=modelMatrix,
                                modelMatrixType=modelMatrixType,
                                quiet=TRUE)
    } else {
      dds.shr <- do.call(rbind, bplapply(levels(parallelIdx), function(l) {
        nbinomWaldTest(dds[parallelIdx == l,,drop=FALSE],
                       betaPrior=TRUE,
                       betaPriorVar=betaPriorVar,
                       modelMatrix=modelMatrix,
                       modelMatrixType=modelMatrixType,
                       quiet=TRUE)
      }, BPPARAM=BPPARAM))
    }
    if (missing(contrast)) {
      # parallel not necessary here
      res.shr <- results(dds.shr, name=coefAlpha, lfcThreshold=lfcThreshold)
    } else {
      # parallel may be useful here as novel contrasts can take a while with big designs
      res.shr <- results(dds.shr, contrast=contrast, lfcThreshold=lfcThreshold,
                         parallel=parallel, BPPARAM=BPPARAM)
    }

    if (lfcThreshold > 0) {
      change.cols <- c("log2FoldChange","lfcSE","stat","pvalue","padj")
    } else {
      change.cols <- c("log2FoldChange","lfcSE")
    }
    for (column in change.cols) {
      res[[column]] <- res.shr[[column]]
    }
    mcols(res,use.names=TRUE)[change.cols,"description"] <- mcols(res.shr,use.names=TRUE)[change.cols,"description"]
    
    deseq2.version <- packageVersion("DESeq2")
    # stash lfcThreshold
    metadata(res)[["lfcThreshold"]] <- lfcThreshold
    priorInfo(res) <- list(type="normal",
                           package="DESeq2",
                           version=deseq2.version,
                           betaPriorVar=betaPriorVar)

    res <- resultsFormatSwitch(object=dds, res=res, format=format)
    return(res)
    
  } else if (type == "apeglm") {

    ############
    ## apeglm ##
    ############
    
    if (!requireNamespace("apeglm", quietly=TRUE)) {
      stop("type='apeglm' requires installing the Bioconductor package 'apeglm'")
    }
    if (!missing(contrast)) {
      stop("type='apeglm' shrinkage only for use with 'coef'")
    }
    stopifnot(!missing(coef))    
    incomingCoef <- gsub(" ","_",sub("log2 fold change \\(MLE\\): ","",mcols(res)[2,2]))
    if (coefAlpha != incomingCoef) {
      stop("'coef' should specify same coefficient as in results 'res'")
    }

    if (!quiet) message("using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    sequence count data: removing the noise and preserving large differences.
    Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895")
    
    Y <- counts(dds)
    if (attr(dds, "modelMatrixType") == "user-supplied") {
      design <- attr(dds, "modelMatrix")
    } else {
      design <- model.matrix(design(dds), data=colData(dds))
    }
    disps <- dispersions(dds)
    if (is.null(normalizationFactors(dds))) {
      offset <- matrix(log(sizeFactors(dds)),
                       nrow=nrow(dds), ncol=ncol(dds), byrow=TRUE)
    } else {
      offset <- log(normalizationFactors(dds))
    }
    if ("weights" %in% assayNames(dds)) {
      weights <- assays(dds)[["weights"]]
    } else {
      weights <- matrix(1, nrow=nrow(dds), ncol=ncol(dds))
    }
    if (apeAdapt) {
      mle <- log(2) * cbind(res$log2FoldChange, res$lfcSE)
    } else {
      mle <- NULL
    }
    if (apeMethod == "general") {
      log.lik <- apeglm::logLikNB
    } else {
      log.lik <- NULL
    }
    if (lfcThreshold > 0) {
      message(paste0("computing FSOS 'false sign or small' s-values (T=",round(lfcThreshold,3),")"))
      svalue <- TRUE
      apeT <- log(2) * lfcThreshold
    } else {
      apeT <- NULL
    }

    # parallel fork
    if (!parallel) {
      fit <- apeglm::apeglm(Y=Y, x=design, log.lik=log.lik, param=disps,
                            coef=coefNum, mle=mle, threshold=apeT,
                            weights=weights, offset=offset,
                            method=apeMethod, ...)
    } else {
      fitList <- bplapply(levels(parallelIdx), function(l) {
        idx <- parallelIdx == l
        apeglm::apeglm(Y=Y[idx,,drop=FALSE], x=design, log.lik=log.lik, param=disps[idx],
                       coef=coefNum, mle=mle, threshold=apeT,
                       weights=weights[idx,,drop=FALSE], offset=offset[idx,,drop=FALSE],
                       method=apeMethod, ...)
      })
      # collate the objects from the split
      fit <- list()
      ape.cols <- c("map","sd","fsr","svalue","interval","diag")
      if (lfcThreshold > 0) {
        ape.cols <- c(ape.cols, "thresh")
      }
      for (param in ape.cols) {
        fit[[param]] <- do.call(rbind, lapply(fitList, `[[`, param))
      }
      fit$prior.control <- fitList[[1]]$prior.control
      fit$svalue <- apeglm::svalue(fit$fsr[,1])
    }
    
    stopifnot(nrow(fit$map) == nrow(dds))
    conv <- fit$diag[,"conv"]
    if (!all(conv[!is.na(conv)] == 0)) {
      message("Some rows did not converge in finding the MAP")
    }
    res$log2FoldChange <- log2(exp(1)) * fit$map[,coefNum]
    res$lfcSE <- log2(exp(1)) * fit$sd[,coefNum]
    mcols(res)$description[2] <- sub("MLE","MAP",mcols(res)$description[2])
    mcols(res)$description[3] <- sub("standard error","posterior SD",mcols(res)$description[3])
    if (svalue) {
      coefAlphaSpaces <- gsub("_"," ",coefAlpha)
      res <- res[,1:3]
      if (lfcThreshold > 0) {
        res$svalue <- as.numeric(apeglm::svalue(fit$thresh))
        description <- paste0("FSOS s-value (T=",round(lfcThreshold,3),"): ",coefAlphaSpaces)
      } else {
        res$svalue <- as.numeric(fit$svalue)
        description <- paste0("s-value: ",coefAlphaSpaces)
      }
      mcols(res)[4,] <- DataFrame(type="results", description=description)
    } else {
      res <- res[,c(1:3,5:6)]
    }
    # stash lfcThreshold
    metadata(res)[["lfcThreshold"]] <- lfcThreshold
    priorInfo(res) <- list(type="apeglm",
                           package="apeglm",
                           version=packageVersion("apeglm"),
                           prior.control=fit$prior.control)
    res <- resultsFormatSwitch(object=dds, res=res, format=format)
    if (returnList) {
      return(list(res=res, fit=fit))
    } else {
      return(res)
    }

  } else if (type == "ashr") {

    ##########
    ## ashr ##
    ##########
    
    if (!requireNamespace("ashr", quietly=TRUE)) {
      stop("type='ashr' requires installing the CRAN package 'ashr'")
    }
    
    if (!quiet) message("using 'ashr' for LFC shrinkage. If used in published research, please cite:
    Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    https://doi.org/10.1093/biostatistics/kxw041")
    
    if (lfcThreshold > 0) message("lfcThreshold is not used by type='ashr'")
    betahat <- res$log2FoldChange
    sebetahat <- res$lfcSE
    fit <- ashr::ash(betahat, sebetahat,
                     mixcompdist="normal", method="shrink", ...)
    res$log2FoldChange <- fit$result$PosteriorMean
    res$lfcSE <- fit$result$PosteriorSD
    mcols(res)$description[2] <- sub("MLE","MMSE",mcols(res)$description[2])
    mcols(res)$description[3] <- sub("standard error","posterior SD",mcols(res)$description[3])
    if (svalue) {
      coefAlphaSpaces <- sub(".*p-value: ","",mcols(res)$description[5])
      res <- res[,1:3]
      res$svalue <- fit$result$svalue
      mcols(res)[4,] <- DataFrame(type="results",
                                  description=paste("s-value:",coefAlphaSpaces))
    } else {
      res <- res[,c(1:3,5:6)]
    }
    # stash lfcThreshold
    metadata(res)[["lfcThreshold"]] <- lfcThreshold
    priorInfo(res) <- list(type="ashr",
                           package="ashr",
                           version=packageVersion("ashr"),
                           fitted_g=fit$fitted_g)
    res <- resultsFormatSwitch(object=dds, res=res, format=format)
    if (returnList) {
      return(list(res=res, fit=fit))
    } else{
      return(res)
    }

  }
}
