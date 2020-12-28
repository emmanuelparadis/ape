\name{dbd}
\alias{dyule}
\alias{dbd}
\alias{dbdTime}
\title{Probability Density Under Birth--Death Models}
\description{
  These functions compute the probability density under some
  birth--death models, that is the probability of obtaining \emph{x}
  species after a time \emph{t} giving how speciation and extinction
  probabilities vary through time (these may be constant, or even equal
  to zero for extinction).
}
\usage{
dyule(x, lambda = 0.1, t = 1, log = FALSE)
dbd(x, lambda, mu, t, conditional = FALSE, log = FALSE)
dbdTime(x, birth, death, t, conditional = FALSE,
        BIRTH = NULL, DEATH = NULL, fast = FALSE)
}
\arguments{
  \item{x}{a numeric vector of species numbers (see Details).}
  \item{lambda}{a numerical value giving the probability of speciation;
    can be a vector with several values for \code{dyule}.}
  \item{mu}{id. for extinction.}
  \item{t}{id. for the time(s).}
  \item{log}{a logical value specifying whether the probabilities should
    be returned log-transformed; the default is \code{FALSE}.}
  \item{conditional}{a logical specifying whether the probabilities
    should be computed conditional under the assumption of no extinction
    after time \code{t}.}
  \item{birth, death}{a (vectorized) function specifying how the
    speciation or extinction probability changes through time (see
    \code{\link{yule.time}} and below).}
  \item{BIRTH, DEATH}{a (vectorized) function giving the primitive
    of \code{birth} or \code{death}.}
  \item{fast}{a logical value specifying whether to use faster
    integration (see \code{\link{bd.time}}).}
}
\details{
  These three functions compute the probabilities to observe \code{x}
  species starting from a single one after time \code{t} (assumed to be
  continuous). The first function is a short-cut for the second one with
  \code{mu = 0} and with default values for the two other arguments.
  \code{dbdTime} is for time-varying \code{lambda} and \code{mu}
  specified as \R functions.

  \code{dyule} is vectorized simultaneously on its three arguments
  \code{x}, \code{lambda}, and \code{t}, according to \R's rules of
  recycling arguments. \code{dbd} is vectorized simultaneously \code{x}
  and \code{t} (to make likelihood calculations easy), and
  \code{dbdTime} is vectorized only on \code{x}; the other arguments are
  eventually shortened with a warning if necessary.

  The returned value is, logically, zero for values of \code{x} out of
  range, i.e., negative or zero for \code{dyule} or if \code{conditional
  = TRUE}. However, it is not checked if the values of \code{x} are
  positive non-integers and the probabilities are computed and returned.

  The details on the form of the arguments \code{birth}, \code{death},
  \code{BIRTH}, \code{DEATH}, and \code{fast} can be found in the links
  below.
}
\note{
  If you use these functions to calculate a likelihood function, it is
  strongly recommended to compute the log-likelihood with, for instance
  in the case of a Yule process, \code{sum(dyule( , log = TRUE))} (see
  examples).
}
\value{
  a numeric vector.
}
\references{
  Kendall, D. G. (1948) On the generalized ``birth-and-death''
  process. \emph{Annals of Mathematical Statistics}, \bold{19}, 1--15.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{bd.time}},  \code{\link{yule.time}}
}
\examples{
x <- 0:10
plot(x, dyule(x), type = "h", main = "Density of the Yule process")
text(7, 0.85, expression(list(lambda == 0.1, t == 1)))

y <- dbd(x, 0.1, 0.05, 10)
z <- dbd(x, 0.1, 0.05, 10, conditional = TRUE)
d <- rbind(y, z)
colnames(d) <- x
barplot(d, beside = TRUE, ylab = "Density", xlab = "Number of species",
        legend = c("unconditional", "conditional on\nno extinction"),
        args.legend = list(bty = "n"))
title("Density of the birth-death process")
text(17, 0.4, expression(list(lambda == 0.1, mu == 0.05, t == 10)))

\dontrun{
### generate 1000 values from a Yule process with lambda = 0.05
x <- replicate(1e3, Ntip(rlineage(0.05, 0)))

### the correct way to calculate the log-likelihood...:
sum(dyule(x, 0.05, 50, log = TRUE))

### ... and the wrong way:
log(prod(dyule(x, 0.05, 50)))

### a third, less preferred, way:
sum(log(dyule(x, 0.05, 50)))
}}
\keyword{utilities}
