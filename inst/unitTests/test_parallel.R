test_parallel <- function() {
  dds <- makeExampleDESeqDataSet()
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersionsGeneEst(dds)
  dds <- estimateDispersionsFit(dds)
  dispPriorVar <- estimateDispersionsPriorVar(dds)
  dds <- estimateDispersionsMAP(dds, dispPriorVar=dispPriorVar)
  dds <- estimateMLEForBetaPriorVar(dds)
  betaPriorVar <- estimateBetaPriorVar(dds)
  dds <- nbinomWaldTest(dds, betaPriorVar=betaPriorVar)
  results(dds)
}
