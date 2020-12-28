\name{corClasses}
\alias{corClasses}
\alias{corPhyl}
\title{Phylogenetic Correlation Structures}
\description{
  Classes of phylogenetic correlation structures (\code{"corPhyl"})
  available in \pkg{ape}.

  \itemize{
    \item{corBrownian}{Brownian motion model (Felsenstein 1985)}
    \item{corMartins}{The covariance matrix defined in Martins and Hansen
      (1997)}
    \item{corGrafen}{The covariance matrix defined in Grafen (1989)}
    \item{corPagel}{The covariance matrix defined in Freckelton et al. (2002)}
    \item{corBlomberg}{The covariance matrix defined in Blomberg et al. (2003)}
  }

  See the help page of each class for references and detailed
  description.
}
\seealso{
  \code{\link[nlme]{corClasses}} and \code{\link[nlme]{gls}} in the
  \pkg{nlme} librarie, \code{\link{corBrownian}},
  \code{\link{corMartins}}, \code{\link{corGrafen}},
  \code{\link{corPagel}}, \code{\link{corBlomberg}},
  \code{\link{vcv}}, \code{\link{vcv2phylo}}
}
\author{Julien Dutheil \email{dutheil@evolbio.mpg.de}, Emmanuel
  Paradis}
\examples{
library(nlme)
txt <- "((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);"
tree.primates <- read.tree(text = txt)
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
Species <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
dat <- data.frame(Species = Species, X = X, Y = Y)

m1 <- gls(Y ~ X, dat, correlation=corBrownian(1, tree.primates, form = ~Species))
summary(m1)
m2 <- gls(Y ~ X, dat, correlation=corMartins(1, tree.primates, form = ~Species))
summary(m2)
corMatrix(m2$modelStruct$corStruct)
m3 <- gls(Y ~ X, dat, correlation=corGrafen(1, tree.primates, form = ~Species))
summary(m3)
corMatrix(m3$modelStruct$corStruct)
}
\keyword{models}

