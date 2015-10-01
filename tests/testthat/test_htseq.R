dir <- system.file(package="pasilla", "extdata")
files <- grep("treated",list.files(dir),value=TRUE)
sampleTable <- data.frame(id=seq_along(files), files,
                          condition=factor(rep(c("t","u"),c(3,4))))
setwd(dir)
expect_error(DESeqDataSetFromHTSeqCount(sampleTable))
dds <- DESeqDataSetFromHTSeqCount(sampleTable, design=~condition)
