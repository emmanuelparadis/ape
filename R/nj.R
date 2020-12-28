## nj.R (2020-09-07)

##   Neighbor-Joining Tree Estimation

## Copyright 2004-2020 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

nj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    if (anyNA(X))
        stop("missing values are not allowed in the distance matrix\nConsider using njs()")
    if (any(is.infinite(X)))
        stop("infinite values are not allowed in the distance matrix")
    N <- as.integer(attr(X, "Size"))
    if (N < 3) stop("cannot build an NJ tree with less than 3 observations")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    obj <- .Call(C_nj, X, N)
    names(obj) <- c("edge", "edge.length")
    dim(obj[[1]]) <- c(2L * N - 3L, 2L)
    obj$tip.label <- labels
    obj$Nnode <- N - 2L
    class(obj) <- "phylo"
    reorder(obj)
}
