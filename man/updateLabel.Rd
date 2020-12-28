\name{updateLabel}
\alias{updateLabel}
\alias{updateLabel.DNAbin}
\alias{updateLabel.AAbin}
\alias{updateLabel.character}
\alias{updateLabel.phylo}
\alias{updateLabel.evonet}
\alias{updateLabel.data.frame}
\alias{updateLabel.matrix}
\title{Update Labels}
\description{
  This function changes labels (names or rownames) giving two vectors (\code{old} and \code{new}). It is a generic function with several methods as described below.
}
\usage{
updateLabel(x, old, new, ...)
\method{updateLabel}{character}(x, old, new, exact = TRUE, ...)
\method{updateLabel}{DNAbin}(x, old, new, exact = TRUE, ...)
\method{updateLabel}{AAbin}(x, old, new, exact = TRUE, ...)
\method{updateLabel}{phylo}(x, old, new, exact = TRUE, nodes = FALSE, ...)
\method{updateLabel}{evonet}(x, old, new, exact = TRUE, nodes = FALSE, ...)
\method{updateLabel}{data.frame}(x, old, new, exact = TRUE, ...)
\method{updateLabel}{matrix}(x, old, new, exact = TRUE, ...)
}
\arguments{
  \item{x}{an object where to change the labels.}
  \item{old, new}{two vectors of mode character (must be of the same length).}
  \item{exact}{a logical value (see details).}
  \item{nodes}{a logical value specifying whether to also update the node labels of the tree or network.}
  \item{\dots}{further arguments passed to and from methods.}
}
\details{
This function can be used to change some of the labels (see examples) or all of them if their ordering is not sure.

If \code{exact = TRUE} (the default), the values in \code{old} are matched exactly with the labels; otherwise (\code{exact = FALSE}), the values in \code{old} are considered as regular expressions and searched in the labels with \code{\link{grep}}.
}
\value{
  an object of the same class than \code{x}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{makeLabel}}, \code{\link{makeNodeLabel}},
  \code{\link{mixedFontLabel}}, \code{\link{stripLabel}},
  \code{\link{checkLabel}}
}
\examples{
\dontrun{
## the tree by Nyakatura & Bininda-Emonds (2012, BMC Biology)
x <- "https://static-content.springer.com/esm/art"
y <- "3A10.1186"
z <- "2F1741-7007-10-12/MediaObjects/12915_2011_534_MOESM5_ESM.NEX"
## The commande below may not print correctly in HTML because of the
## percentage symbol; see the text or PDF help page.
url <- paste(x, y, z, sep = "%")
TC <- read.nexus(url)
tr <- TC$carnivoreST_bestEstimate
old <- c("Uncia_uncia", "Felis_manul", "Leopardus_jacobitus")
new <- c("Panthera_uncia", "Otocolobus_manul", "Leopardus_jacobita")
tr.updated <- updateLabel(tr, old, new)
}

tr <- rtree(6)
## the order of the labels are randomized by this function
old <- paste0("t", 1:6)
new <- paste0("x", 1:6)
updateLabel(tr, old, new)
tr
}
\keyword{manip}
