## me.R (2019-03-26)

##   Tree Estimation Based on Minimum Evolution Algorithm

## Copyright 2007 Vincent Lefort with modifications by
##                Emmanuel Paradis (2008-2019)

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

fastme.bal <- function(X, nni = TRUE, spr = TRUE, tbr = FALSE)
{
    if (tbr) {
        warning("option 'tbr = TRUE' was ignored: see ?fastme.bal")
        tbr <- FALSE
    }
    if (is.matrix(X)) X <- as.dist(X)
    N <- as.integer(attr(X, "Size"))
    nedge <- 2L * N - 3L
    ans <- .C(me_b, as.double(X), N, 1:N, as.integer(nni),
              as.integer(spr), as.integer(tbr), integer(nedge),
              integer(nedge), double(nedge), NAOK = TRUE)
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    labels <- labels[ans[[3]]]
    obj <- list(edge =  cbind(ans[[7]], ans[[8]]),
                edge.length = ans[[9]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"
    obj
}

fastme.ols <- function(X, nni = TRUE)
{
    if (is.matrix(X)) X <- as.dist(X)
    N <- as.integer(attr(X, "Size"))
    nedge <- 2L * N - 3L
    ans <- .C(me_o, as.double(X), N, 1:N, as.integer(nni),
              integer(nedge), integer(nedge), double(nedge),
              NAOK = TRUE)
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    labels <- labels[ans[[3]]]
    obj <- list(edge =  cbind(ans[[5]], ans[[6]]),
                edge.length = ans[[7]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    attr(obj, "order") <- "cladewise"
    obj
}

bionj <- function(X)
{
    if (is.matrix(X)) X <- as.dist(X)
    if (any(is.na(X)))
        stop("missing values are not allowed in the distance matrix.\nConsider using bionjs()")
    if (any(X > 100))
        stop("at least one distance was greater than 100")
    N <- as.integer(attr(X, "Size"))

    ans <- .C(C_bionj, as.double(X), N, integer(2 * N - 3),
              integer(2 * N - 3), double(2*N - 3), NAOK = TRUE)
    labels <- attr(X, "Labels")
    if (is.null(labels)) labels <- as.character(1:N)
    obj <- list(edge =  cbind(ans[[3]], ans[[4]]), edge.length = ans[[5]],
                tip.label = labels, Nnode = N - 2L)
    class(obj) <- "phylo"
    reorder(obj)
}
