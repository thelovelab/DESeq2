makeExpandedModelMatrix <- function(object) {
  designFactors <- getDesignFactors(object)
  coldata <- colData(object)
  coldata <- rbind(coldata,coldata[nrow(coldata),])
  for (f in designFactors) {
    levels(coldata[[f]]) <- c(levels(coldata[[f]]),"_null_level_")
    coldata[[f]] <- relevel(coldata[[f]],"_null_level_")
    coldata[[f]][nrow(coldata)] <- "_null_level_"
  }
  mm0 <- model.matrix(design(object), data=coldata)
  mm <- mm0[-nrow(mm0),]
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  colnames(mm) <- make.names(colnames(mm))
  mm
}

averagePriorsOverLevels <- function(object, betaPriorVar) { 
  expandedModelMatrix <- makeExpandedModelMatrix(object)
  expandedNames <- colnames(expandedModelMatrix)
  betaPriorIn <- betaPriorVar
  betaPriorOut <- numeric(length(expandedNames))
  bpiNms <- names(betaPriorIn)
  idx <- which(bpiNms %in% expandedNames)
  betaPriorOut[match(bpiNms[idx],expandedNames)] <- betaPriorIn[idx]
  designFactors <- getDesignFactors(object)
  coldata <- colData(object)
  for (f in designFactors) {
    lvls <- levels(coldata[[f]])
    mmColnames <- make.names(paste0(f,c(lvls,"Cntrst")))
    meanPriorVar <- mean(betaPriorIn[names(betaPriorIn) %in% mmColnames])
    betaPriorOut[expandedNames %in% mmColnames] <- meanPriorVar
  }
  # also set prior for any interactions between design factors
  # which are new in the expanded model matrix using existing interactions
  termsOrder <- attr(terms.formula(design(object)),"order")
  if (any(termsOrder > 1)) {
    for (f1 in designFactors) {
      for (f2 in designFactors) {
        if (f1 == f2) next
        lvls1 <- levels(coldata[[f1]])
        lvls2 <- levels(coldata[[f2]])
        mmColnames <- make.names(paste0(f1,rep(lvls1,each=length(lvls2)),":",
                                        f2,rep(lvls2,times=length(lvls1))))
        meanPriorVar <- mean(betaPriorIn[names(betaPriorIn) %in% mmColnames])
        betaPriorOut[expandedNames %in% mmColnames] <- meanPriorVar
      }
    }
  }
  betaPriorOut
}

# adds all first order contrasts
addAllContrasts <- function(object, betaMatrix) { 
  designFactors <- getDesignFactors(object)
  coldata <- colData(object)
  for (f in designFactors) {
    lvls <- levels(coldata[[f]])
    mmColnames <- make.names(paste0(f,lvls))
    M <- betaMatrix[,colnames(betaMatrix) %in% mmColnames,drop=FALSE]
    n <- ncol(M)
    if (n > 1) {
      if (n == 2) {
        is <- 2
        js <- 1
      } else {
        is <- do.call(c,sapply(seq_len(n-1)+1, function(k) seq(from=k,to=n)))
        js <- rep(seq_len(n-1),rev(seq_len(n-1)))
      }
      contrastCols <- mapply(function(i,j) M[,i] - M[,j], i=is, j=js)
      colnames(contrastCols) <- rep(make.names(paste0(f,"Cntrst")),ncol(contrastCols))
      betaMatrix <- cbind(betaMatrix, contrastCols)
    }
  }
  betaMatrix
}
  
# want to make a model matrix which won't change from releveling
# this function is only for use in calculating beta prior
# in the case of expanded model matrices, which shouldn't change with releveling
makeReleveledModelMatrix <- function(object) {
  designFactors <- getDesignFactors(object)
  coldata <- colData(object)
  # pick an arbitrary sample for setting base levels
  # either the sample with smallest size factor, or the first sample
  sf <- sizeFactors(object)
  if (!is.null(sf)) idx <- which.min(sf) else idx <- 1
  for (v in designFactors) {
    coldata[[v]] <- relevel(coldata[[v]], as.character(coldata[[v]][idx]))
  }
  mm <- model.matrix(design(object), data=coldata)
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  colnames(mm) <- make.names(colnames(mm))
  mm
}

