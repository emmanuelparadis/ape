## treePop.R (2011-10-11)

## Tree Reconstruction With the Triangles Method

## Copyright 2011 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

triangMtd <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    if (any(is.na(X)))
        stop("missing values are not allowed in the distance matrix")
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_triangMtd, as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3), NAOK = TRUE)
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}

triangMtds <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_triangMtds, as.double(X), as.integer(N), integer(2*N - 3),
              integer(2*N - 3), double(2*N - 3), NAOK = TRUE)
    obj <- list(edge = cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
