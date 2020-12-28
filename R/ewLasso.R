## ewLasso.R (2013-04-02)

##   Lasso Tree

## Copyright 2013 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ewLasso <- function(X, phy)
{
    if (is.matrix(X)) X <- as.dist(X)
    X[is.na(X)] <- -1
    X[X < 0] <- -1
    X[is.nan(X)] <- -1

    if (is.rooted(phy)) {
        phy <- unroot(phy)
        warning("'phy' is rooted: it was unrooted for this operation")
    }

    N <- attr(X, "Size")
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    ans <- .C(C_ewLasso, as.double(X), as.integer(N),
              as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
              NAOK = TRUE)
}
