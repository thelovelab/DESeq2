\name{results}
\alias{removeResults}
\alias{results}
\alias{resultsNames}
\title{Extract results from a DESeq analysis}
\usage{
  results(object, name)

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
}
\value{
  For \code{results}: a DataFrame of results columns with
  metadata columns of coefficient and test information

  For \code{resultsNames}: the names of the columns
  available as results, usually a combination of the
  variable name and a level

  For \code{removeResults}: the original object with
  results metadata columns removed
}
\description{
  Extract results from a DESeq analysis giving base means
  across samples, log2 fold changes, standard errors,
  p-values and adjusted p-values. Find available names for
  results; Return an object with results columns removed.
}
\details{
  Multiple results can be returned for an analysis with
  multifactor design, so \code{results} takes an argument
  \code{name}

  The available names can be checked using
  \code{resultsNames}; these are combined variable names
  and factor levels, potentially with minor changes made by
  the \code{DataFrame} function on column names (e.g.
  dashes into periods).

  By default, results for the last variable will be
  returned. Information on the variable represented and the
  test used for p-values (Wald test or likelihood ratio
  test) is stored in the metadata columns, accessible by
  calling \code{mcol} on the object returned by
  \code{results}.

  For analyses using the likelihood ratio test (using
  \code{\link{nbinomLRT}}), the p-values are determined
  solely by the difference in deviance between the full and
  reduced model formula. In this case, the \code{name}
  argument only specifies which coefficient should be used
  for reporting the log2 fold changes.

  Results can be removed from an object by calling
  \code{removeResults}
}
\examples{
example("DESeq")
results(dds)
resultsNames(dds)
dds <- removeResults(dds)
}
\seealso{
  \code{\link{DESeq}}
}

