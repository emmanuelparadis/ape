## unique.multiPhylo.R (2021-12-19)

##   Revomes Duplicate Trees from a List

## Copyright 2007-2021 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

unique.multiPhylo <-
    function(x, incomparables = FALSE,
             use.edge.length = FALSE,
             use.tip.label = TRUE, ...)
{
    n <- length(x)
    ## fixed by Martin:
    if (n == 0L) return(x)
    if (n == 1L) return(structure(x, old.index = 1L))
    keep <- 1L
    
    tmp <- try(.compressTipLabel(x), TRUE)
    if(requireNamespace("phangorn") && !inherits(tmp, "try-error") && !use.edge.length && use.tip.label){
        hash_x <- phangorn::hash(x)
        keep <- which(!duplicated(hash_x))
        old.index <- match(hash_x, unique(hash_x))
    } else {
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
    }
    res <- x[keep]
    attr(res, "old.index") <- old.index
    res
}
