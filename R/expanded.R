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
  # these can appear when interactions are present without main effect variables
  nullLvls <- grepl("_null_level_",colnames(mm0))
  mm <- mm0[-nrow(mm0), !nullLvls]
  attr(mm,"assign") <- attr(mm0,"assign")
  colnames(mm)[colnames(mm) == "(Intercept)"] <- "Intercept"
  colnames(mm) <- make.names(colnames(mm))
  mm
}

averagePriorsOverLevels <- function(object, betaPriorVar) { 
  expandedModelMatrix <- makeExpandedModelMatrix(object)
  expandedNames <- colnames(expandedModelMatrix)
  betaPriorIn <- betaPriorVar
  betaPriorOut <- numeric(length(expandedNames))
  names(betaPriorOut) <- expandedNames
  bpiNms <- names(betaPriorIn)
  idx <- which(bpiNms %in% expandedNames)
  betaPriorOut[match(bpiNms[idx],expandedNames)] <- betaPriorIn[idx]
  designFactors <- getDesignFactors(object)
  allVars <- all.vars(design(object))
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
      for (f2 in allVars) {
        if (f1 == f2) next
        lvls1 <- levels(coldata[[f1]])
        # the case where f2 is a factor like f1
        if (f2 %in% designFactors) {
          lvls2 <- levels(coldata[[f2]])
          mmColnames <- make.names(paste0(f1,rep(lvls1,each=length(lvls2)),":",
                                          f2,rep(lvls2,times=length(lvls1))))
          meanPriorVar <- mean(betaPriorIn[names(betaPriorIn) %in% mmColnames])
          betaPriorOut[expandedNames %in% mmColnames] <- meanPriorVar
        # the case where f2 is not a factor
        } else {
          mmColnames <- make.names(c(paste0(f1,lvls1,":",f2),paste0(f2,":",f1,lvls1)))
          meanPriorVar <- mean(betaPriorIn[names(betaPriorIn) %in% mmColnames])
          betaPriorOut[expandedNames %in% mmColnames] <- meanPriorVar
        }
      }
    }
  }
  if (any(is.na(betaPriorOut))) {
    stop(paste("beta prior for",paste(names(betaPriorOut)[is.na(betaPriorOut)],collapse=","),"is NA"))
  }
  if (!all(betaPriorOut > 0)) {
    stop(paste("beta prior for",paste(names(betaPriorOut)[betaPriorOut <= 0],collapse=","),"is not greater than 0"))
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
  

