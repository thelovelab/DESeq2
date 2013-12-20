test_LRT_prior <- function() {
  dds <- makeExampleDESeqDataSet(n=100)
  colData(dds)$condition <- factor(rep(1:3,each=4))
  colData(dds)$group <- factor(rep(1:2,6))
  design(dds) <- ~ group + condition
  dds <- DESeq(dds,test="LRT",reduced=~ group,betaPrior=TRUE)
  checkTrue(any(attr(dds,"betaPriorVar") < 1e6))
  res <- results(dds)
  checkTrue(grepl("LRT",mcols(res)$description[colnames(res) == "stat"]))

  design(dds) <- ~ group * condition
  dds <- DESeq(dds,test="LRT",reduced=~ group,betaPrior=TRUE)
  res <- results(dds)
}
