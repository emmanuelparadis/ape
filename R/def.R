## def.R (2014-10-24)

##   Definition of Vectors for Plotting or Annotating

## Copyright 2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

def <- function(x, ..., default = NULL, regexp = FALSE)
{
    dots <- list(...)
    if (is.null(default)) {
        if (is.numeric(dots[[1L]])) default <- 1
        if (is.character(dots[[1L]])) default <- "black"
    }
    foo <- if (regexp) function(vec, y) grep(y, vec) else function(vec, y) which(vec == y)
    res <- rep(default, length(x))
    nms <- names(dots)
    for (i in seq_along(nms))
        res[foo(x, nms[i])] <- dots[[i]]
    res
}
