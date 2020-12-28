## mvr.R (2012-03-30)

##   Minimum Variance Reduction

## Copyright 2011 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

mvr <- function(X, V)
{
    if (is.matrix(X)) X <- as.dist(X)
    if (is.matrix(V)) V <- as.dist(V)
    if (any(is.na(X)))
        stop("missing values are not allowed in the distance matrix")
    if (any(is.na(V)))
        stop("missing values are not allowed in the variance matrix")
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_mvr, as.double(X), as.double(V), as.integer(N),
              integer(2*N - 3), integer(2*N - 3), double(2*N - 3),
              NAOK = TRUE)
    obj <- list(edge = cbind(ans[[4]], ans[[5]]), edge.length = ans[[6]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}

mvrs <- function(X, V, fs = 15)
{
    if (fs < 1)
        stop("argument 'fs' must be a non-zero positive integer")

    if (is.matrix(X)) X <- as.dist(X)
    if (is.matrix(V)) V <- as.dist(V)

    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1

    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_mvrs, as.double(X), as.double(V), as.integer(N),
              integer(2*N - 3), integer(2*N - 3), double(2*N - 3),
              as.integer(fs), NAOK = TRUE)
    obj <- list(edge = cbind(ans[[4]], ans[[5]]), edge.length = ans[[6]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
