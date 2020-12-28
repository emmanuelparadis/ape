\name{compar.cheverud}
\alias{compar.cheverud}
\title{Cheverud's Comparative Method}
\description{
  This function computes the phylogenetic variance component and the
  residual deviation for continous characters, taking into account the
  phylogenetic relationships among species, following the comparative
  method described in Cheverud et al. (1985). The correction proposed by
  Rholf (2001) is used.
}
\usage{
compar.cheverud(y, W, tolerance = 1e-06, gold.tol = 1e-04)
}
\arguments{
  \item{y}{A vector containing the data to analyse.}
  \item{W}{The phylogenetic connectivity matrix. All diagonal elements
    will be ignored.}
  \item{tolerance}{Minimum difference allowed to consider eigenvalues as
    distinct.}
  \item{gold.tol}{Precision to use in golden section search alogrithm.}
}
\details{
  Model: \deqn{y = \rho W y + e}{y = rho.W.y + e}

  where \eqn{e}{e} is the error term, assumed to be normally distributed.
  \eqn{\rho}{rho} is estimated by the maximum likelihood procedure given
  in Rohlf (2001), using a golden section search algorithm. The code of
  this function is indeed adapted from a MatLab code given in appendix
  in Rohlf's article, to correct a mistake in Cheverud's original paper.
}
\value{
  A list with the following components:

  \item{rhohat}{The maximum likelihood estimate of \eqn{\rho}{rho}}
  \item{Wnorm}{The normalized version of \code{W}}
  \item{residuals}{Error terms (\eqn{e}{e})}
}
\references{
  Cheverud, J. M., Dow, M. M. and Leutenegger, W. (1985) The quantitative
  assessment of phylogenetic constraints in comparative analyses: sexual
  dimorphism in body weight among primates. \emph{Evolution},
  \bold{39}, 1335--1351.

  Rohlf, F. J. (2001) Comparative methods for the analysis of continuous
  variables: geometric interpretations. \emph{Evolution}, \bold{55},
  2143--2160.

  Harvey, P. H. and Pagel, M. D. (1991) \emph{The Comparative Method in
    Evolutionary Biology}. Oxford University Press.
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de}}
\seealso{\code{\link{compar.lynch}}}
\examples{
### Example from Harvey and Pagel's book:
y<-c(10,8,3,4)
W <- matrix(c(1,1/6,1/6,1/6,1/6,1,1/2,1/2,1/6,1/2,1,1,1/6,1/2,1,1), 4)
compar.cheverud(y,W)
### Example from Rohlf's 2001 article:
W<- matrix(c(
  0,1,1,2,0,0,0,0,
  1,0,1,2,0,0,0,0,
  1,1,0,2,0,0,0,0,
  2,2,2,0,0,0,0,0,
  0,0,0,0,0,1,1,2,
  0,0,0,0,1,0,1,2,
  0,0,0,0,1,1,0,2,
  0,0,0,0,2,2,2,0
),8)
W <- 1/W
W[W == Inf] <- 0
y<-c(-0.12,0.36,-0.1,0.04,-0.15,0.29,-0.11,-0.06)
compar.cheverud(y,W)
}
\keyword{regression}
