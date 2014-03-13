\docType{methods}
\name{fpkm}
\alias{fpkm}
\title{FPKM: fragments per kilobase per million mapped fragments}
\usage{
fpkm(object, robust = TRUE)
}
\arguments{
  \item{object}{a \code{DESeqDataSet}}

  \item{robust}{whether to use size factors to normalize
  rather than taking the column sums of the raw counts,
  using the \code{\link{fpm}} function.}
}
\value{
a matrix which is normalized per kilobase of the union of
basepairs in the \code{GRangesList} or \code{GRanges} of
the mcols(object), and per million of mapped fragments,
either using the robust median ratio method (robust=TRUE,
default) or using raw counts (robust=FALSE). If the
\code{mcol(object)} is not a \code{GRangesList} or
\code{GRanges}, a column
\code{mcols(object)$unionBasepairs} should be defined for
the kilobase normalization.
}
\description{
FPKM: fragments per kilobase per million mapped fragments
}
\examples{
# create a matrix with 1 million counts for the
# 2nd and 3rd column, the 1st and 4th have
# half and double the counts, respectively.
m <- matrix(1e6 * rep(c(.125, .25, .25, .5), each=4),
            ncol=4, dimnames=list(1:4,1:4))
mode(m) <- "integer"
se <- SummarizedExperiment(m, colData=DataFrame(sample=1:4))
dds <- DESeqDataSet(se, ~ 1)

# create 4 GRanges with lengths: 1, 1, 2, 2.5 Kb
gr1 <- GRanges("chr1",IRanges(1,1000))
gr2 <- GRanges("chr1",IRanges(c(1,501),c(500,1000)))
gr3 <- GRanges("chr1",IRanges(c(1,1001),c(1000,2000)))
gr4 <- GRanges("chr1",IRanges(c(1,1001,2001),c(500,3000,3000)))
rowData(dds) <- GRangesList(gr1,gr2,gr3,gr4)

# the raw counts
counts(dds)

# the FPKM values
fpkm(dds)

# held constant per 1 million fragments
counts(dds) <- counts(dds) * 2L
round(fpkm(dds))
}
\seealso{
\code{\link{fpm}}
}

