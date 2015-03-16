dds <- makeExampleDESeqDataSet(n=100, m=18)
dds$group <- factor(rep(1:3,each=6))
dds$condition <- factor(rep(rep(c("A","B","C"),each=2),3))
# note: design is not used
design(dds) <- ~ 1
dds <- dds[,-c(17,18)]

m1 <- model.matrix(~ group*condition, colData(dds))
m1 <- m1[,-9]
m0 <- model.matrix(~ group + condition, colData(dds))

dds <- DESeq(dds, full=m1, reduced=m0, test="LRT")
results(dds)[1,]
results(dds, name="group2.conditionC", test="Wald")[1,]
dds <- removeResults(dds)
dds <- DESeq(dds, full=m1, test="Wald", betaPrior=FALSE)
results(dds)[1,]

