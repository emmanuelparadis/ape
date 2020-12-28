\name{def}
\alias{def}
\title{Definition of Vectors for Plotting or Annotating}
\description{
  This function can be used to define vectors to annotate a set of taxon
  names, labels, etc. It should facilitate the (re)definition of colours
  or similar attributes for plotting trees or other graphics.
}
\usage{
def(x, ..., default = NULL, regexp = FALSE)
}
\arguments{
  \item{x}{a vector of mode character.}
  \item{\dots}{a series of statements defining the attributes.}
  \item{default}{the default to be used (see details).}
  \item{regexp}{a logical value specifying whether the statements
    defined in \code{\dots} should be taken as regular expressions.}
}
\details{
The idea of this function is to make the definition of colours, etc.,
simpler than what is done usually. A typical use is:

\code{def(tr$tip.label, Homo_sapiens = "blue")}

which will return a vector of character strings all "black" except one
matching the tip label "Homo_sapiens" which will be "blue". Another use
could be:

\code{def(tr$tip.label, Homo_sapiens = 2)}

which will return a vector a numerical values all 1 except for
"Homo_sapiens" which will be 2. Several definitions can be done, e.g.:

\code{def(tr$tip.label, Homo_sapiens = "blue", Pan_paniscus = "red")}

The default value is determined with respect to the mode of the values
given with the \code{\dots} (either "black" or 1).

If \code{regexp = TRUE} is used, then the names of the statements must be
quoted, e.g.:

\code{def(tr$tip.label, "^Pan_" = "red", regexp = TRUE)}

will return "red" for all labels starting with "Pan_".
}
\value{
a vector of the same length than \code{x}.
}
\author{Emmanuel Paradis}
\examples{
data(bird.orders)
a <- def(bird.orders$tip.label, Galliformes = 2)
str(a) # numeric
plot(bird.orders, font = a)
co <- def(bird.orders$tip.label, Passeriformes = "red", Trogoniformes = "blue")
str(co) # character
plot(bird.orders, tip.color = co)
### use of a regexp (so we need to quote it) to colour all orders
### with names starting with "C" (and change the default):
co2 <- def(bird.orders$tip.label, "^C" = "gold", default = "grey", regexp = TRUE)
plot(bird.orders, tip.color = co2)
}
\keyword{manip}
