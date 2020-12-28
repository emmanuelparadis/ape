\name{diversity.contrast.test}
\alias{diversity.contrast.test}
\title{Diversity Contrast Test}
\description{
  This function performs the diversity contrast test comparing pairs of
  sister-clades.
}
\usage{
diversity.contrast.test(x, method = "ratiolog",
        alternative = "two.sided", nrep = 0, ...)
}
\arguments{
  \item{x}{a matrix or a data frame with at least two columns: the first
    one gives the number of species in clades with a trait supposed to
    increase or decrease diversification rate, and the second one the
    number of species in the sister-clades without the trait. Each row
    represents a pair of sister-clades.}
  \item{method}{a character string specifying the kind of test:
    \code{"ratiolog"} (default), \code{"proportion"},
    \code{"difference"}, \code{"logratio"}, or any unambiguous
    abbreviation of these.}
  \item{alternative}{a character string defining the alternative
    hypothesis:  \code{"two.sided"} (default), \code{"less"},
    \code{"greater"}, or any unambiguous abbreviation of these.}
  \item{nrep}{the number of replications of the randomization test; by
    default, a Wilcoxon test is done.}
  \item{\dots}{arguments passed to the function
    \code{\link[stats]{wilcox.test}}.}
}
\details{
  If \code{method = "ratiolog"}, the test described in Barraclough et
  al. (1996) is performed. If \code{method = "proportion"}, the version
  in Barraclough et al. (1995) is used. If \code{method = "difference"},
  the signed difference is used (Sargent 2004). If \code{method = "logratio"},
  then this is Wiegmann et al.'s (1993) version. These
  four tests are essentially different versions of the same test (Vamosi
  and Vamosi 2005, Vamosi 2007). See Paradis (2012) for a comparison of
  their statistical performance with other tests.

  If \code{nrep = 0}, a Wilcoxon test is done on the species diversity
  contrasts with the null hypothesis is that they are distributed around
  zero. If \code{nrep > 0}, a randomization procedure is done where the
  signs of the diversity contrasts are randomly chosen. This is used to
  create a distribution of the test statistic which is compared with the
  observed value (the sum of the diversity contrasts).
}
\value{
  a single numeric value with the \emph{P}-value.
}
\references{
  Barraclough, T. G., Harvey, P. H. and  Nee, S. (1995) Sexual
  selection and taxonomic diversity in passerine birds.
  \emph{Proceedings of the Royal Society of London. Series B. Biological
    Sciences}, \bold{259}, 211--215.

  Barraclough, T. G., Harvey, P. H., and Nee,  S. (1996) Rate of
  \emph{rbc}L gene sequence evolution and species diversification in
  flowering plants (angiosperms). \emph{Proceedings of the Royal Society
    of London. Series B. Biological Sciences}, \bold{263}, 589--591.

  Paradis, E. (2012) Shift in diversification in sister-clade
  comparisons: a more powerful test. \emph{Evolution}, \bold{66},
  288--295.

  Sargent, R. D. (2004) Floral symmetry affects speciation rates in
  angiosperms. \emph{Proceedings of the Royal Society of London. Series
  B. Biological Sciences}, \bold{271}, 603--608.

  Vamosi, S. M. (2007) Endless tests: guidelines for analysing non-nested
  sister-group comparisons. An addendum. \emph{Evolutionary Ecology
  Research}, \bold{9}, 717.

  Vamosi, S. M. and Vamosi, J. C. (2005) Endless tests: guidelines for
  analysing non-nested sister-group comparisons. \emph{Evolutionary
  Ecology Research}, \bold{7}, 567--579.

  Wiegmann, B., Mitter, C. and Farrell, B. 1993. Diversification of
  carnivorous parasitic insects: extraordinary radiation or specialized
  dead end? \emph{American Naturalist}, \bold{142}, 737--754.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{slowinskiguyer.test}}, \code{\link{mcconwaysims.test}}
  \code{\link{richness.yule.test}}
}
\examples{
### data from Vamosi & Vamosi (2005):
fleshy <- c(1, 1, 1, 1, 1, 3, 3, 5, 9, 16, 33, 40, 50, 100, 216, 393, 850, 947,1700)
dry <- c(2, 64, 300, 89, 67, 4, 34, 10, 150, 35, 2, 60, 81, 1, 3, 1, 11, 1, 18)
x <- cbind(fleshy, dry)
diversity.contrast.test(x)
diversity.contrast.test(x, alt = "g")
diversity.contrast.test(x, alt = "g", nrep = 1e4)
slowinskiguyer.test(x)
mcconwaysims.test(x)
}
\keyword{htest}
