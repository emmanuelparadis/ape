\name{all.equal.DNAbin}
\alias{all.equal.DNAbin}
\title{Compare DNA Sets}
\description{
  Comparison of DNA sequence sets, particularly when aligned.
}
\usage{
\method{all.equal}{DNAbin}(target, current, plot = FALSE, ...)
}
\arguments{
  \item{target, current}{the two sets of sequences to be compared.}
  \item{plot}{a logical value specifying whether to plot the sites that
    are different (only if the labels of both alignments are the same).}
  \item{\dots}{further arguments passed to \code{\link{image.DNAbin}}.}
}
\details{
  If the two sets of DNA sequences are exactly identical, this function
  returns \code{TRUE}. Otherwise, a detailed comparison is made only if
  the labels (i.e., rownames) of \code{target} and \code{current} are the
  same (possibly in different orders). In all other cases, a brief
  description of the differences is returned (sometimes with
  recommendations to make further comparisons).

  This function can be used for testing in programs using
  \code{\link[base]{isTRUE}} (see examples below).
}
\value{
  \code{TRUE} if the two sets are identical; a list with two elements
  (message and different.sites) if a detailed comparison is done; or a
  vector of mode character.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{image.DNAbin}}, \code{\link{clustal}},
  \code{\link{checkAlignment}},
  the generic function: \code{\link[base]{all.equal}}
}
\examples{
data(woodmouse)
woodm2 <- woodmouse
woodm2[1, c(1:5, 10:12, 30:40)] <- as.DNAbin("g")
res <- all.equal(woodmouse, woodm2, plot = TRUE)
str(res)

## if used for testing in R programs:
isTRUE(all.equal(woodmouse, woodmouse)) # TRUE
isTRUE(all.equal(woodmouse, woodm2)) # FALSE

all.equal(woodmouse, woodmouse[15:1, ])
all.equal(woodmouse, woodmouse[-1, ])
all.equal(woodmouse, woodmouse[, -1])

\dontrun{
## To run the followings you need internet and Clustal and MUSCLE
## correctly installed.
## Data from Johnson et al. (2006, Science)
refs <- paste("DQ082", 505:545, sep = "")
DNA <- read.GenBank(refs)
DNA.clustal <- clustal(DNA)
DNA.muscle <- muscle(DNA)
isTRUE(all.equal(DNA.clustal, DNA.muscle)) # FALSE
all.equal(DNA.clustal, DNA.muscle, TRUE)
}
}
\keyword{manip}
