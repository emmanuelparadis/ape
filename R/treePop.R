## treePop.R (2011-10-11)

##   Tree Popping

## Copyright 2011 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

treePop <- function(obj)
{
    mf <- obj$matsplit
    labels <- obj$labels
    n <- length(labels)
    imf <- as.integer(mf)
    freq <- obj$freq
    mimf <- matrix(imf, nrow(mf), ncol(mf))
    ans <- .C(C_treePop, mimf, as.double(freq), as.integer(ncol(mf)),
              as.integer(n), integer(2*n - 3), integer(2*n - 3),
              double(2*n - 3), NAOK = TRUE)
    obj <- list(edge = cbind(ans[[5]], ans[[6]]), edge.length = ans[[7]],
                tip.label = labels, Nnode = n - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
