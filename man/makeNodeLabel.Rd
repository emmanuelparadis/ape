\name{makeNodeLabel}
\alias{makeNodeLabel}
\title{Makes Node Labels}
\description{
  This function makes node labels in a tree in a flexible way.
}

\usage{
makeNodeLabel(phy, method = "number", prefix = "Node", nodeList = list(), ...)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{method}{a character string giving the method used to create the
    labels. Three choices are possible: \code{"number"} (the default),
    \code{"md5sum"}, and \code{"user"}, or any unambiguous abbreviation
    of these.}
  \item{prefix}{the prefix used if \code{method = "number"}.}
  \item{nodeList}{a named list specifying how nodes are names if
    \code{method = "user"} (see details and examples).}
   \item{\dots}{further arguments passed to \code{grep}.}
}
\details{
  The three methods are described below:

  \itemize{
    \item{``number''}{The labels are created with 1, 2, \dots prefixed
      with the argument \code{prefix}; thus the default is to have
      Node1, Node2, \dots Set \code{prefix = ""} to have only numbers.}
    \item{``md5sum''}{For each node, the labels of the tips descendant
      from this node are extracted, sorted alphabetically, and written
      into a temporary file, then the md5sum of this file is extracted
      and used as label. This results in a 32-character string which is
      unique (even accross trees) for a given set of tip labels.}
    \item{``user''}{the argument \code{nodeList} must be a list with
      names, the latter will be used as node labels. For each element of
      \code{nodeList}, the tip labels of the tree are searched for
      patterns present in this element: this is done using
      \code{\link[base]{grep}}. Then the most recent common ancestor of
      the matching tips is given the corresponding names as labels. This
      is repeated for each element of \code{nodeList}.}
  }

  The method \code{"user"} can be used in combination with either of the
  two others (see examples). Note that this method only modifies the
  specified node labels (so that if the other nodes have already labels
  they are not modified) while the two others change all labels.
}
\value{
  an object of class \code{"phylo"}.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{makeLabel}}, \code{\link[base]{grep}},
  \code{\link{mixedFontLabel}}, \code{\link{label2table}},
  \code{\link{checkLabel}}
}
\examples{
tr <-
"((Pan_paniscus,Pan_troglodytes),((Homo_sapiens,Homo_erectus),Homo_abilis));"
tr <- read.tree(text = tr)
tr <- makeNodeLabel(tr, "u", nodeList = list(Pan = "Pan", Homo = "Homo"))
plot(tr, show.node.label = TRUE)
### does not erase the previous node labels:
tr <- makeNodeLabel(tr, "u", nodeList = list(Hominid = c("Pan","Homo")))
plot(tr, show.node.label = TRUE)
### the two previous commands could be combined:
L <- list(Pan = "Pan", Homo = "Homo", Hominid = c("Pan","Homo"))
tr <- makeNodeLabel(tr, "u", nodeList = L)
### combining different methods:
tr <- makeNodeLabel(tr, c("n", "u"), prefix = "#", nodeList = list(Hominid = c("Pan","Homo")))
plot(tr, show.node.label = TRUE)
}
\keyword{manip}
