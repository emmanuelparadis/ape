\name{compar.lynch}
\alias{compar.lynch}
\title{Lynch's Comparative Method}
\usage{
compar.lynch(x, G, eps = 1e-4)
}
\arguments{
  \item{x}{eiher a matrix, a vector, or a data.frame containing the data
    with species as rows and variables as columns.}
  \item{G}{a matrix that can be interpreted as an among-species correlation
    matrix.}
  \item{eps}{a numeric value to detect convergence of the EM algorithm.}
}
\description{
  This function computes the heritable additive value and the residual
  deviation for continous characters, taking into account the
  phylogenetic relationships among species, following the comparative
  method described in Lynch (1991).
}
\details{
  The parameter estimates are computed following the EM
  (expectation-maximization) algorithm. This algorithm usually leads to
  convergence but may lead to local optima of the likelihood
  function. It is recommended to run several times the function in order
  to detect these potential local optima. The `optimal' value for
  \code{eps} depends actually on the range of the data and may be
  changed by the user in order to check the stability of the parameter
  estimates. Convergence occurs when the differences between two
  successive iterations of the EM algorithm leads to differences between
  both residual and additive values less than or equal to \code{eps}.
}
\note{
  The present function does not perform the estimation of ancestral
  phentoypes as proposed by Lynch (1991). This will be implemented in
  a future version.
}
\value{
  A list with the following components:
  \item{vare}{estimated residual variance-covariance matrix.}
  \item{vara}{estimated additive effect variance covariance matrix.}
  \item{u}{estimates of the phylogeny-wide means.}
  \item{A}{addtitive value estimates.}
  \item{E}{residual values estimates.}
  \item{lik}{logarithm of the likelihood for the entire set of observed
    taxon-specific mean.}
}
\references{
  Lynch, M. (1991) Methods for the analysis of comparative data in
  evolutionary biology. \emph{Evolution}, \bold{45}, 1065--1080.
}
\author{Julien Claude \email{julien.claude@umontpellier.fr}}
\seealso{
  \code{\link{pic}}, \code{\link{compar.gee}}
}
\examples{
### The example in Lynch (1991)
cat("((((Homo:0.21,Pongo:0.21):0.28,",
   "Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);",
   file = "ex.tre", sep = "\n")
tree.primates <- read.tree("ex.tre")
unlink("ex.tre")
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
compar.lynch(cbind(X, Y),
             G = vcv.phylo(tree.primates, cor = TRUE))
}
\keyword{regression}
