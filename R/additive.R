## additive.R (2013-10-04)

##   Incomplete Distance Matrix Filling

## Copyright 2011-2013 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

additive <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    m <- sum(X == -1)
    ans <- .C(C_additive, as.double(X), as.integer(N),
              as.integer(m), double(N*N))
    matrix(ans[[4]], N, N)
}

ultrametric <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1
    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    m <- sum(X == -1)
    ans <- .C(C_ultrametric, as.double(X), as.integer(N),
              as.integer(m), double(N*N))
    matrix(ans[[4]], N, N)
}
