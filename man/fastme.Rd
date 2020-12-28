\name{FastME}
\alias{FastME}
\alias{fastme}
\alias{fastme.bal}
\alias{fastme.ols}
\title{
  Tree Estimation Based on the Minimum Evolution Algorithm
}
\description{
  The two FastME functions (balanced and OLS) perform the
  minimum evolution algorithm of Desper and Gascuel (2002).
}
\usage{
  fastme.bal(X, nni = TRUE, spr = TRUE, tbr = FALSE)
  fastme.ols(X, nni = TRUE)
}
\arguments{
  \item{X}{a distance matrix; may be an object of class \code{"dist"}.}
  \item{nni}{a logical value; TRUE to perform NNIs (default).}
  \item{spr}{ditto for SPRs.}
  \item{tbr}{ignored (see details).}
}
\details{
  The code to perform topology searches based on TBR (tree bisection and
  reconnection) did not run correctly and has been removed after the
  release of \pkg{ape} 5.3. A warning is issued if \code{tbr = TRUE}.
}
\value{
  an object of class \code{"phylo"}.
}
\references{
  Desper, R. and Gascuel, O. (2002) Fast and accurate phylogeny
  reconstruction algorithms based on the minimum-evolution principle.
  \emph{Journal of Computational Biology}, \bold{9}, 687--705.
}
\author{
  original C code by Richard Desper; adapted and ported to R
  by Vincent Lefort \email{vincent.lefort@lirmm.fr}
}
\seealso{
  \code{\link{nj}}, \code{\link{bionj}},
  \code{\link{write.tree}}, \code{\link{read.tree}},
  \code{\link{dist.dna}}
}
\examples{
### From Saitou and Nei (1987, Table 1):
x <- c(7, 8, 11, 13, 16, 13, 17, 5, 8, 10, 13,
       10, 14, 5, 7, 10, 7, 11, 8, 11, 8, 12,
       5, 6, 10, 9, 13, 8)
M <- matrix(0, 8, 8)
M[lower.tri(M)] <- x
M <- t(M)
M[lower.tri(M)] <- x
dimnames(M) <- list(1:8, 1:8)
tr <- fastme.bal(M)
plot(tr, "u")
### a less theoretical example
data(woodmouse)
trw <- fastme.bal(dist.dna(woodmouse))
plot(trw)
}
\keyword{models}

