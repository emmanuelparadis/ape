## unique.multiPhylo.R (2026-01-19)

##   Revomes Duplicate Trees from a List

## Copyright 2007-2021 Emmanuel Paradis, 2025 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

hash_topo <- function(x, ...){
  x <- unroot(x)
  x <- reorder(x, "postorder")
  fun <- function(x){
    nTips <- as.integer(length(x$tip.label))
    res <- sorted_bipartition(x$edge, nTips)
    l <- lengths(res)
    res <- res[l>1]
    digest(res, ...)
  }
  if(inherits(x, "phylo")) return(fun(x, ...))
  if(inherits(x, "multiPhylo")){
    x <- .compressTipLabel(x)
    return(sapply(x, fun, ...))
  }
  stop("x needs to be an object of class phylo or multiPhylo!")
}


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
    if(!inherits(tmp, "try-error") && !use.edge.length && use.tip.label){
        hash_x <- hash_topo(x)
        keep <- which(!duplicated(hash_x))
        old.index <- match(hash_x, hash_x[keep])
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
