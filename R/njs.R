## njs.R (2013-10-04)

## Tree Reconstruction from Incomplete Distances With NJ* or bio-NJ*

## Copyright 2011-2013 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

njs <- function(X, fs = 15)
{
    if (fs < 1)
        stop("argument 'fs' must be a non-zero positive integer")
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_njs, as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3), as.integer(fs),
              NAOK = TRUE)
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}

bionjs <- function(X, fs = 15)
{
    if (fs < 1)
        stop("argument 'fs' must be a non-zero positive integer")
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_bionjs, as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3), as.integer(fs),
              NAOK = TRUE)
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
