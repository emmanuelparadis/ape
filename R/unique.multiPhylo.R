## unique.multiPhylo.R (2014-01-15)

##   Revomes Duplicate Trees from a List

## Copyright 2007-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

unique.multiPhylo <-
    function(x, incomparables = FALSE,
             use.edge.length = FALSE,
             use.tip.label = TRUE, ...)
{
    n <- length(x)
    keep <- 1L
    old.index <- seq_len(n)
    for (i in 2:n) {
        already.seen <- FALSE
        for (j in keep) {
            if (all.equal(x[[j]], x[[i]],
                          use.edge.length = use.edge.length,
                          use.tip.label = use.tip.label)) {
                already.seen <- TRUE
                old.index[i] <- j
                break
            }
        }
        if (!already.seen) keep <- c(keep, i)
    }
    res <- x[keep]
    attr(res, "old.index") <- old.index
    res
}
