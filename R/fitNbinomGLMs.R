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
                          useQR=TRUE, forceOptim=FALSE, warnNonposVar=TRUE, minmu=0.5,
                          type = c("DESeq2", "glmGamPoi")) {
  type <- match.arg(type, c("DESeq2", "glmGamPoi"))
  
  if (missing(modelFormula)) {
    modelFormula <- design(object)
  }
  if (is.null(modelMatrix)) {
    modelAsFormula <- TRUE
    modelMatrix <- stats::model.matrix.default(modelFormula, data=as.data.frame(colData(object)))
  } else {
    modelAsFormula <- FALSE
  }

  stopifnot(all(MatrixGenerics::colSums(abs(modelMatrix)) > 0))

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
  
  normalizationFactors <- getSizeOrNormFactors(object)
  
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

  # use weights if they are present in assays(object)
  wlist <- getAndCheckWeights(object, modelMatrix)
  weights <- wlist$weights
  useWeights <- wlist$useWeights
  
  if(type == "glmGamPoi"){
    stopifnot("type = 'glmGamPoi' cannot handle weights" = ! useWeights,
              "type = 'glmGamPoi' does not support NA's in alpha_hat" = all(! is.na(alpha_hat))) 
    gp_res <- glmGamPoi::glm_gp(counts(object), design = modelMatrix,
                                size_factors = FALSE, offset = log(normalizationFactors),
                                overdispersion = alpha_hat, verbose = FALSE)
    logLikeMat <- dnbinom(counts(object), mu=gp_res$Mu, size=1/alpha_hat, log=TRUE)
    logLike <- MatrixGenerics::rowSums(logLikeMat)
    res <- list(logLike = logLike, betaConv =  rep(TRUE, nrow(object)), betaMatrix = gp_res$Beta / log(2),
                betaSE = NULL, mu = gp_res$Mu, betaIter = rep(NA,nrow(object)),
                modelMatrix=modelMatrix, 
                nterms=ncol(modelMatrix), hat_diagonals = NULL)
    return(res)
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
      betaMatrix <- if (useWeights) {
                      matrix(log2(MatrixGenerics::rowSums(weights*counts(object, normalized=TRUE))
                                  /MatrixGenerics::rowSums(weights)),ncol=1)
                    } else {
                      matrix(log2(MatrixGenerics::rowMeans(counts(object, normalized=TRUE))),ncol=1)
                    }
      mu <- normalizationFactors * as.numeric(2^betaMatrix)
      logLikeMat <- dnbinom(counts(object), mu=mu, size=1/alpha, log=TRUE)
      logLike <- if (useWeights) {
                   MatrixGenerics::rowSums(weights*logLikeMat)
                 } else {
                   MatrixGenerics::rowSums(logLikeMat)
                 }
      modelMatrix <- stats::model.matrix.default(~ 1, data=as.data.frame(colData(object)))
      colnames(modelMatrix) <- modelMatrixNames <- "Intercept"
      w <- if (useWeights) {
             weights * (mu^-1 + alpha)^-1
           } else {
             (mu^-1 + alpha)^-1
           }
      xtwx <- MatrixGenerics::rowSums(w)
      sigma <- xtwx^-1
      betaSE <- matrix(log2(exp(1)) * sqrt(sigma),ncol=1)      
      hat_diagonals <- w * xtwx^-1;
      res <- list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
                  betaSE = betaSE, mu = mu, betaIter = betaIter,
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
      logBaseMean <- log(MatrixGenerics::rowMeans(counts(object,normalized=TRUE)))
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
  lambdaNatLogScale <- lambda / log(2)^2
  
  betaRes <- fitBetaWrapper(ySEXP = counts(object), xSEXP = modelMatrix,
                            nfSEXP = normalizationFactors,
                            alpha_hatSEXP = alpha_hat,
                            beta_matSEXP = beta_mat,
                            lambdaSEXP = lambdaNatLogScale,
                            weightsSEXP = weights,
                            useWeightsSEXP = useWeights,
                            tolSEXP = betaTol, maxitSEXP = maxit,
                            useQRSEXP=useQR, minmuSEXP=minmu)

  # Note on deviance: the 'deviance' calculated in fitBeta() (C++)
  # is not returned in mcols(object)$deviance. instead, we calculate
  # the log likelihood below and use -2 * logLike.
  # (reason is that we have other ways of estimating beta:
  # above intercept code, and below optim code)
  
  mu <- normalizationFactors * t(exp(modelMatrix %*% t(betaRes$beta_mat)))
  dispersionVector <- rep(dispersions(object), times=ncol(object))
  logLike <- nbinomLogLike(counts(object), mu, dispersions(object), weights, useWeights)

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
                                   weights,useWeights,
                                   betaMatrix,betaSE,betaConv,
                                   beta_mat,
                                   mu,logLike,minmu=minmu)
    betaMatrix <- resOptim$betaMatrix
    betaSE <- resOptim$betaSE
    betaConv <- resOptim$betaConv
    mu <- resOptim$mu
    logLike <- resOptim$logLike
  }

  stopifnot(!any(is.na(betaSE)))
  nNonposVar <- sum(MatrixGenerics::rowSums(betaSE == 0) > 0)
  if (warnNonposVar & nNonposVar > 0) warning(nNonposVar,"rows had non-positive estimates of variance for coefficients")
  
  list(logLike = logLike, betaConv = betaConv, betaMatrix = betaMatrix,
       betaSE = betaSE, mu = mu, betaIter = betaRes$iter, modelMatrix=modelMatrix, 
       nterms=ncol(modelMatrix), hat_diagonals=betaRes$hat_diagonals)
}

# this function calls fitNbinomGLMs() twice:
# 1 - without the beta prior, in order to calculate the
#     beta prior variance and hat matrix
# 2 - again but with the prior in order to get beta matrix and standard errors
fitGLMsWithPrior <- function(object, betaTol, maxit, useOptim, useQR, betaPriorVar, modelMatrix=NULL, minmu=0.5) {
  
  objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
  modelMatrixType <- attr(object, "modelMatrixType")

  if (missing(betaPriorVar) | !(all(c("mu","H") %in% assayNames(objectNZ)))) {

    # stop unless modelMatrix was NOT supplied, the code below all works
    # by building model matrices using the formula, doesn't work with incoming model matrices
    stopifnot(is.null(modelMatrix))
    
    # fit the negative binomial GLM without a prior,
    # used to construct the prior variances
    # and for the hat matrix diagonals for calculating Cook's distance
    fit <- fitNbinomGLMs(objectNZ,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         renameCols = (modelMatrixType == "standard"),
                         minmu=minmu)
    modelMatrix <- fit$modelMatrix
    modelMatrixNames <- colnames(modelMatrix)
    H <- fit$hat_diagonals
    betaMatrix <- fit$betaMatrix
    mu <- fit$mu

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
    if (is.null(modelMatrix)) {
      modelMatrix <- getModelMatrix(object)
    }
    H <- assays(objectNZ)[["H"]]
    mu <- assays(objectNZ)[["mu"]]
    mleBetaMatrix <- as.matrix(mcols(objectNZ)[,grep("MLE_",names(mcols(objectNZ))),drop=FALSE])
  }
     
  if (missing(betaPriorVar)) {
    betaPriorVar <- estimateBetaPriorVar(objectNZ, modelMatrix=modelMatrix)
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
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         minmu=minmu)
    modelMatrix <- fit$modelMatrix
  } else if (modelMatrixType == "expanded") {
    modelMatrix <- makeExpandedModelMatrix(objectNZ)
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         modelMatrix=modelMatrix, renameCols=FALSE,
                         minmu=minmu)
  } else if (modelMatrixType == "user-supplied") {
    fit <- fitNbinomGLMs(objectNZ, lambda=lambda,
                         betaTol=betaTol, maxit=maxit,
                         useOptim=useOptim, useQR=useQR,
                         modelMatrix=modelMatrix, renameCols=FALSE,
                         minmu=minmu)
  }

  res <- list(fit=fit, H=H, betaPriorVar=betaPriorVar, mu=mu,
              modelMatrix=modelMatrix, mleBetaMatrix=mleBetaMatrix)
  res
}

# breaking out the optim backup code from fitNbinomGLMs
fitNbinomGLMsOptim <- function(object,modelMatrix,lambda,
                               rowsForOptim,rowStable,
                               normalizationFactors,alpha_hat,
                               weights,useWeights,
                               betaMatrix,betaSE,betaConv,
                               beta_mat,
                               mu,logLike,minmu=0.5) {
  x <- modelMatrix
  lambdaNatLogScale <- lambda / log(2)^2
  large <- 30
  for (row in rowsForOptim) {
    betaRow <- if (rowStable[row] & all(abs(betaMatrix[row,]) < large)) {
      betaMatrix[row,]
    } else {
      beta_mat[row,]
    }
    nf <- normalizationFactors[row,]
    k <- counts(object)[row,]
    alpha <- alpha_hat[row]
    objectiveFn <- function(p) {
      mu_row <- as.numeric(nf * 2^(x %*% p))
      logLikeVector <- dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE)
      logLike <- if (useWeights) {
                   sum(weights[row,] * logLikeVector)
                 } else {
                   sum(logLikeVector)
                 }
      logPrior <- sum(dnorm(p,0,sqrt(1/lambda),log=TRUE))
      negLogPost <- -1 * (logLike + logPrior)
      if (is.finite(negLogPost)) negLogPost else 10^300
    }
    o <- optim(betaRow, objectiveFn, method="L-BFGS-B",lower=-large, upper=large)
    ridge <- if (length(lambdaNatLogScale) > 1) {
      diag(lambdaNatLogScale)
    } else {
      as.matrix(lambdaNatLogScale,ncol=1)
    }
    # if we converged, change betaConv to TRUE
    if (o$convergence == 0) {
      betaConv[row] <- TRUE
    }
    # with or without convergence, store the estimate from optim
    betaMatrix[row,] <- o$par
    # calculate the standard errors
    mu_row <- as.numeric(nf * 2^(x %*% o$par))
    # store the new mu vector
    mu[row,] <- mu_row
    mu_row[mu_row < minmu] <- minmu
    w <- if (useWeights) {
           diag(weights[row,] * (mu_row^-1 + alpha)^-1)
         } else {
           diag((mu_row^-1 + alpha)^-1)
         }
    xtwx <- t(x) %*% w %*% x
    xtwxRidgeInv <- solve(xtwx + ridge)
    sigma <- xtwxRidgeInv %*% xtwx %*% xtwxRidgeInv
    # warn below regarding these rows with negative variance
    betaSE[row,] <- log2(exp(1)) * sqrt(pmax(diag(sigma),0))
    logLikeVector <- dnbinom(k,mu=mu_row,size=1/alpha,log=TRUE)
    logLike[row] <- if (useWeights) {
                      sum(weights[row,] * logLikeVector)
                    } else {
                      sum(logLikeVector)
                    }
  }
  return(list(betaMatrix=betaMatrix,betaSE=betaSE,
              betaConv=betaConv,mu=mu,logLike=logLike))
}
