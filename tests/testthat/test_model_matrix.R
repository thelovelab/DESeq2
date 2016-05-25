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

# test better error than "error: inv(): matrix seems singular"
coldata <- data.frame(group=factor(rep(1:3,each=6)),
                      group2=factor(rep(1:3,each=6)),
                      condition=factor(rep(1:6,3)))
counts <- matrix(rpois(180, 100), ncol=18)
m1 <- model.matrix(~ group + group2, coldata)
m2 <- model.matrix(~ condition + group, coldata)
dds <- DESeqDataSetFromMatrix(counts, coldata, ~group)
expect_error(dds <- DESeq(dds, full=m1, fitType="mean"), "full rank")
expect_error(dds <- DESeq(dds, full=m2, reduced=m1, test="LRT", fitType="mean"), "full rank")
