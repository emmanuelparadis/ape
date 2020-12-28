\name{matexpo}
\alias{matexpo}
\title{Matrix Exponential}
\usage{
matexpo(x)
}
\arguments{
  \item{x}{a square matrix of mode numeric.}
}
\description{
  This function computes the exponential of a square matrix using a
  spectral decomposition.
}
\value{
  a numeric matrix of the same dimensions than `x'.
}
\author{Emmanuel Paradis}
\examples{
### a simple rate matrix:
m <- matrix(0.1, 4, 4)
diag(m) <- -0.3
### towards equilibrium:
for (t in c(1, 5, 10, 50)) print(matexpo(m*t))
}
\keyword{array}
\keyword{multivariate}
