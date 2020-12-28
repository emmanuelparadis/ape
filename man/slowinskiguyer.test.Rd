\name{slowinskiguyer.test}
\alias{slowinskiguyer.test}
\title{Slowinski-Guyer Test of Homogeneous Diversification}
\description{
  This function performs the Slowinski--Guyer test that a trait or
  variable does not increase diversification rate.
}
\usage{
slowinskiguyer.test(x, detail = FALSE)
}
\arguments{
  \item{x}{a matrix or a data frame with at least two columns: the first
    one gives the number of species in clades with a trait supposed to
    increase diversification rate, and the second one the number of
    species in the corresponding sister-clade without the trait. Each
    row represents a pair of sister-clades.}
  \item{detail}{if \code{TRUE}, the individual P-values are appended.}
}
\details{
  The Slowinski--Guyer test compares a series of sister-clades where one
  of the two is characterized by a trait supposed to increase
  diversification rate. The null hypothesis is that the trait does not
  affect diversification. If the trait decreased diversification rate,
  then the null hypothesis cannot be rejected.

  The present function has mainly a historical interest. The
  Slowinski--Guyer test generally performs poorly: see Paradis (2012)
  alternatives and the functions cited below.
}
\value{
  a data frame with the \eqn{\chi^2}{chi2}, the number of degrees of
  freedom, and the \emph{P}-value. If \code{detail = TRUE}, a list is
  returned with the data frame and a vector of individual
  \emph{P}-values for each pair of sister-clades.
}
\references{
  Paradis, E. (2012) Shift in diversification in sister-clade
  comparisons: a more powerful test. \emph{Evolution}, \bold{66},
  288--295.

  Slowinski, J. B. and Guyer, C. (1993) Testing whether certain traits
  have caused amplified diversification: an improved method based on a
  model of random speciation and extinction. \emph{American Naturalist},
  \bold{142}, 1019--1024.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{balance}}, \code{\link{mcconwaysims.test}},
  \code{\link{diversity.contrast.test}},
  \code{\link{richness.yule.test}},
  \code{\link[geiger]{rc}} in \pkg{geiger},
  \code{\link[apTreeshape]{shift.test}} in \pkg{apTreeshape}
}
\examples{
### from Table 1 in Slowinski and Guyer(1993):
viviparous <- c(98, 8, 193, 36, 7, 128, 2, 3, 23, 70)
oviparous <- c(234, 17, 100, 4, 1, 12, 6, 1, 481, 11)
x <- data.frame(viviparous, oviparous)
slowinskiguyer.test(x, TRUE) # 'P ~ 0.32' in the paper
xalt <- x
xalt[3, 2] <- 1
slowinskiguyer.test(xalt)
}
\keyword{htest}
