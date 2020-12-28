## mst.R (2020-07-02)

##   Minimum Spanning Tree

## Copyright 2002-2006 Yvonnick Noel, Julien Claude, and Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

mst <- function(X)
{
    if (inherits(X, "dist")) X <- as.matrix(X)
    n <- dim(X)[1]
    N <- matrix(0, n, n)
    tree <- NULL
    large.value <- max(X) + 1
    diag(X) <- large.value
    index.i <- 1

    for (i in 1:(n - 1)) {
        tree <- c(tree, index.i)
        m <- apply(as.matrix(X[, tree]), 2, min)  #calcul les minimum par colonne
        a <- sortIndex(X[, tree])[1, ]
        b <- sortIndex(m)[1]
        index.j <- tree[b]
        index.i <- a[b]

        N[index.i, index.j] <- 1
        N[index.j, index.i] <- 1

        for (j in tree) {
            X[index.i, j] <- large.value
            X[j, index.i] <- large.value
        }
    }
    dimnames(N) <- dimnames(X)
    class(N) <- "mst"
    return(N)
}

### Function returning an index matrix for an increasing sort
sortIndex <- function(X)
{
    if(length(X) == 1) return(1)                  # sorting a scalar?
    if(!is.matrix(X)) X <- as.matrix(X)           # force vector into matrix
    ## n <- nrow(X)
    apply(X, 2, function(v) order(rank(v)))       # find the permutation
}

plot.mst <- function(x, graph = "circle", x1 = NULL, x2 = NULL, ...)
{
    n <- nrow(x)
    if (is.null(x1) || is.null(x2)) {
        if (graph == "circle") {
            ang <- seq(0, 2 * pi, length = n + 1)
            x1 <- cos(ang)
            x2 <- sin(ang)
            plot(x1, x2, type = "n", xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n", bty = "n", ...)
        }
        if (graph == "nsca") {
            XY <- nsca(x)
            x1 <- XY[, 1]
            x2 <- XY[, 2]
            plot(XY, type = "n", xlab = "\"nsca\" -- axis 1",
                 ylab = "\"nsca\" -- axis 2", ...)
        }
    } else plot(x1, x2, type = "n", xlab = deparse(substitute(x1)),
                ylab = deparse(substitute(x2)), ...)
    for (i in 1:n) {
        w1 <- which(x[i, ] == 1)
        segments(x1[i], x2[i], x1[w1], x2[w1])
    }
    points(x1, x2, pch = 21, col = "black", bg = "white", cex = 3)
    text(x1, x2, 1:n, cex = 0.8)
}

nsca <- function(A)
{
    Dr <- apply(A, 1, sum)
    Dc <- apply(A, 2, sum)

    eig.res <- eigen(diag(1 / sqrt(Dr)) %*% A %*% diag(1 / sqrt(Dc)))
    r <- diag(1 / Dr) %*% (eig.res$vectors)[, 2:4]
    ## The next line has been changed by EP (20-02-2003), since
    ## it does not work if 'r' has no dimnames already defined
    ## dimnames(r)[[1]] <- dimnames(A)[[1]]
    rownames(r) <- rownames(A)
    r
}
