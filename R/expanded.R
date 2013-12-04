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

averagePriorsOverLevels <- function(object, betaPriorVar, modelMatrixNames) { 
  expandedModelMatrix <- makeExpandedModelMatrix(object)
  expandedNames <- colnames(expandedModelMatrix)
  betaPriorIn <- betaPriorVar
  betaPriorOut <- numeric(length(expandedNames))
  betaPriorOut[match(modelMatrixNames,expandedNames)] <- betaPriorIn
  designFactors <- getDesignFactors(object)
  coldata <- colData(object)
  for (f in designFactors) {
    lvls <- levels(coldata[[f]])
    mmColnames <- make.names(paste0(f,lvls))
    meanPriorVar <- mean(betaPriorIn[modelMatrixNames %in% mmColnames])
    betaPriorOut[expandedNames %in% mmColnames] <- meanPriorVar
  }
  betaPriorOut
}

    
