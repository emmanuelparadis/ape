## ladderize.R (2007-10-08)

##   Matrix Exponential

## Copyright 2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

matexpo <- function(x)
{
    if (!is.matrix(x)) stop('"x" must be a matrix')
    nr <- dim(x)[1]
    if (nr != dim(x)[2]) stop('"x" must be a square matrix')
    ans <- .C(mat_expo, as.double(x), as.integer(nr))[[1]]
    dim(ans) <- c(nr, nr)
    ans
}
