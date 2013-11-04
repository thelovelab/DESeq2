\name{results}
\alias{removeResults}
\alias{results}
\alias{resultsNames}
\title{Extract results from a DESeq analysis}
\usage{
  results(object, name, contrast, lfcThreshold = 0,
    altHypothesis = c("greaterAbs", "lessAbs", "greater", "less"),
    cooksCutoff, independentFiltering = TRUE, alpha = 0.1,
    filter, theta = seq(0, 0.95, by = 0.05),
    pAdjustMethod = "BH")

  resultsNames(object)

  removeResults(object)
}
\arguments{
  \item{object}{a DESeqDataSet, on which one of the
  following functions has already been called:
  \code{\link{DESeq}}, \code{\link{nbinomWaldTest}}, or
  \code{\link{nbinomLRT}}}

  \item{name}{the name of the coefficient for which to
  report log2 fold changes -- and for the Wald test,
  p-values and adjusted p-values}

  \item{contrast}{either a character vector with exactly
  three elements (name of the factor, name of the numerator
  level, name of the denominator level), or a numeric
  contrast vector with one element for each element in
  \code{resultsNames(object)}, i.e. columns of the model
  matrix. The DESeqDataSet must be one produced using the
  Wald test steps in order to use contrasts. If used, the
  \code{name} argument is ignored.}

  \item{lfcThreshold}{a non-negative value, which specifies
  the test which should be applied to the log2 fold
  changes. The standard is a test that the log2 fold
  changes are not equal to zero. However, log2 fold changes
  greater or less than \code{lfcThreshold} can also be
  tested. Specify the alternative hypothesis using the
  \code{altHypothesis} argument. If \code{lfcThreshold} is
  specified, the results are Wald tests, and LRT p-values
  will be overwritten.}

  \item{altHypothesis}{character which specifies the
  alternative hypothesis, i.e. those values of log2 fold
  change which the user is interested in finding. The
  complement of this set of values is the null hypothesis
  which will be tested. If the log2 fold change specified
  by \code{name} or by \code{contrast} is written as \eqn{
  \beta }{ beta }, then the possible values for
  \code{altHypothesis} represent the following alternate
  hypotheses: \itemize{ \item greaterAbs - \eqn{|\beta| >
  \textrm{lfcThreshold} }{ |beta| > lfcThreshold }, and
  p-values are two-tailed \item lessAbs - \eqn{ |\beta| <
  \textrm{lfcThreshold} }{ |beta| < lfcThreshold }, and
  p-values are the maximum of the upper and lower tests.
  \item greater - \eqn{ \beta > \textrm{lfcThreshold} }{
  beta > lfcThreshold } \item less - \eqn{ \beta <
  -\textrm{lfcThreshold} }{ beta < -lfcThreshold } } Note:
  for \code{altHypothesis="lessAbs"}, a value less than
  1/(1 + baseMean) for base mean less than 5 is added to
  the standard error, to stabilize the ratio lfc/lfcSE when
  numerator and denominator tend to zero.}

  \item{cooksCutoff}{theshold on Cook's distance, such that
  if one or more samples for a row have a distance higher,
  the p-value for the row is set to NA. The default cutoff
  is the .75 quantile of the F(p, m-p) distribution, where
  p is the number of coefficients being fitted and m is the
  number of samples. Set to Inf or FALSE to disable the
  resetting of p-values to NA. Note: this test excludes the
  Cook's distance of samples whose removal would result in
  rank deficient design matrix.}

  \item{independentFiltering}{logical, whether independent
  filtering should be applied automatically}

  \item{alpha}{the significance cutoff used for optimizing
  the independent filtering}

  \item{filter}{the vector of filter statistics over which
  the independent filtering will be optimized. By default
  the mean of normalized counts is used.}

  \item{theta}{the quantiles at which to assess the number
  of rejections from independent filtering}

  \item{pAdjustMethod}{the method to use for adjusting
  p-values, see \code{?p.adjust}}
}
\value{
  For \code{results}: a DataFrame of results columns:
  \code{log2FoldChange}, \code{lfcSE}, \code{stat},
  \code{pvalue} and \code{padj}. also includes metadata
  columns of variable information

  For \code{resultsNames}: the names of the columns
  available as results, usually a combination of the
  variable name and a level

  For \code{removeResults}: the original object with
  results metadata columns removed
}
\description{
  \code{results} extracts results from a DESeq analysis
  giving base means across samples, log2 fold changes,
  standard errors, test statistics, p-values and adjusted
  p-values. \code{resultsNames} finds available names for
  results; \code{removeResults} returns a
  \code{DESeqDataSet} object with results columns removed.
}
\details{
  Multiple results can be returned for an analysis with
  multifactor design, so \code{results} takes an argument
  \code{name} as well as the argument \code{contrast}
  explained below.

  The available names can be checked using
  \code{resultsNames}; these are combined variable names
  and factor levels, potentially with minor changes made by
  the \code{make.names} function on column names (e.g.
  dashes into periods).

  By default, results for the last variable will be
  returned. Information on the variable represented and the
  test used for p-values (Wald test or likelihood ratio
  test) is stored in the metadata columns, accessible by
  calling \code{mcol} on the \code{DataFrame} returned by
  \code{results}.

  By default, independent filtering is performed to select
  a set of genes which will result in the most genes with
  adjusted p-values less than a threshold, alpha. The
  adjusted p-values for the genes which do not pass the
  filter threshold are set to \code{NA}. By default, the
  mean of normalized counts is used to perform this
  filtering, though other statistics can be provided.
  Several arguments from the \code{filtered_p} function of
  genefilter are provided to control or turn off the
  independent filtering behavior.

  A contrast can be performed by specifying the variable
  and factor levels which should be compared, or by
  specifying the numeric contrast vector. The test
  statistic is then:

  \deqn{ c^t \beta / \sqrt{c^t \Sigma c } }{ c' beta /
  sqrt( c' Sigma c ) }

  For analyses using the likelihood ratio test (using
  \code{\link{nbinomLRT}}), the p-values are determined
  solely by the difference in deviance between the full and
  reduced model formula. In this case, the \code{name}
  argument only specifies which coefficient should be used
  for reporting the log2 fold changes.

  Cook's distance for each sample are accessible as a
  matrix "cooks" stored in the assays() list. This measure
  is useful for identifying rows where the observed counts
  might not fit to a negative binomial distribution.

  Results can be removed from an object by calling
  \code{removeResults}
}
\examples{
example("DESeq")
results(dds)
resultsNames(dds)
dds <- removeResults(dds)
}
\references{
  Richard Bourgon, Robert Gentleman, Wolfgang Huber:
  Independent filtering increases detection power for
  high-throughput experiments. PNAS (2010),
  \url{http://dx.doi.org/10.1073/pnas.0914005107}
}
\seealso{
  \code{\link{DESeq}}
}

