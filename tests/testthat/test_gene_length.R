if (FALSE) {

n <- 10
files <- c("sample1","sample2","sample3","sample4")

envir <- environment()

# Cufflinks-like
gene_id <- rep(paste0("gene",seq_len(n)),each=3)
set.seed(1)
sample1 <- data.frame(gene_id=gene_id,length=rpois(3*n,2000),FPKM=round(rnorm(3*n,10,1),2))
sample2 <- data.frame(gene_id=gene_id,length=rpois(3*n,2000),FPKM=round(rnorm(3*n,10,1),2))
sample3 <- data.frame(gene_id=gene_id,length=rpois(3*n,2000),FPKM=round(rnorm(3*n,10,1),2))
sample4 <- data.frame(gene_id=gene_id,length=rpois(3*n,2000),FPKM=round(rnorm(3*n,10,1),2))
importer <- get
dds <- makeExampleDESeqDataSet(n=n + 2, m=4)
dds <- normalizeGeneLength(dds, files=files, level="tx",
                           geneIdCol="gene_id", lengthCol="length", abundanceCol="FPKM",
                           dropGenes=TRUE, importer=importer, envir=envir)

# RSEM-like
gene_id <- paste0("gene",seq_len(n))
sample1 <- data.frame(gene_id=gene_id,effective_length=rpois(n,2000))
sample2 <- data.frame(gene_id=gene_id,effective_length=rpois(n,2000))
sample3 <- data.frame(gene_id=gene_id,effective_length=rpois(n,2000))
sample4 <- data.frame(gene_id=gene_id,effective_length=rpois(n,2000))
dds <- makeExampleDESeqDataSet(n=n + 2, m=4)
dds <- normalizeGeneLength(dds, files=files, level="gene",
                           geneIdCol="gene_id", lengthCol="effective_length",
                           dropGenes=TRUE, importer=importer, envir=envir)


}
