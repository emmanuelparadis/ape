\name{correlogram.formula}
\alias{correlogram.formula}
\title{Phylogenetic Correlogram}
\usage{
  correlogram.formula(formula, data = NULL, use = "all.obs")
}
\arguments{
  \item{formula}{a formula of the type \code{y1+..+yn ~ g1/../gn}, where
    the \code{y}'s are the data to analyse and the \code{g}'s are the
    taxonomic levels.}
  \item{data}{a data frame containing the variables specified in the
    formula. If \code{NULL}, the variables are sought in the user's
    workspace.}
  \item{use}{a character string specifying how to handle missing
    values (i.e., \code{NA}). This must be one of  "all.obs",
    "complete.obs", or "pairwise.complete.obs", or any unambiguous
    abbrevation of these. In the first case, the presence of missing
    values produces an error. In the second case, all rows with missing
    values will be removed before computation. In the last case, missing
    values are removed on a case-by-case basis.}
}
\description{
  This function computes a correlogram from taxonomic levels.
}
\details{
  See the vignette in R: \code{vignette("MoranI")}.
}
\value{
  An object of class \code{correlogram} which is a data frame with three
  columns:

  \item{obs}{the computed Moran's I}
  \item{p.values}{the corresponding P-values}
  \item{labels}{the names of each level}

  or an object of class \code{correlogramList} containing a list of
  objects of class \code{correlogram} if several variables are given as
  response in \code{formula}.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de} and
  Emmanuel Paradis}
\seealso{
  \code{\link{plot.correlogram}, \link{Moran.I}}
}
\examples{
data(carnivora)
### Using the formula interface:
co <- correlogram.formula(SW ~ Order/SuperFamily/Family/Genus,
      data=carnivora)
co
plot(co)
### Several correlograms on the same plot:
cos <- correlogram.formula(SW + FW ~ Order/SuperFamily/Family/Genus,
      data=carnivora)
cos
plot(cos)
}
\keyword{regression}
