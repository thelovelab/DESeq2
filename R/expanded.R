getDesignFactors <- function(object) {
  design <- design(object)
  designVars <- all.vars(formula(design))
  designVarsClass <- sapply(designVars, function(v) class(colData(object)[[v]]))
  designVars[designVarsClass == "factor"]
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
  mm
}

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
  betaPriorOut
}

addAllContrasts <- function(object, betaMatrix) { 
  designFactors <- getDesignFactors(object)
  coldata <- colData(object)
  for (f in designFactors) {
    lvls <- levels(coldata[[f]])
    mmColnames <- make.names(paste0(f,lvls))
    M <- betaMatrix[,colnames(betaMatrix) %in% mmColnames]
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

    
