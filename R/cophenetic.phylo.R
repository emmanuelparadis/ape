## cophenetic.phylo.R (2024-12-10)

##   Pairwise Distances from a Phylogenetic Tree

## Copyright 2006-2024 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.nodes <- function(x, fail.if.no.length = FALSE)
{
    if (is.null(x$edge.length)) {
        if (fail.if.no.length) stop("the tree has no branch length")
        warning("the tree has no branch length: fixing them to one.")
        x <- compute.brlen(x, 1)
    }
    x <- reorder(x) # required for the C code
    n <- Ntip(x)
    m <- x$Nnode
    d <- .Call(dist_nodes, n, m, x$edge, x$edge.length)
    nm <- n + m
    dimnames(d) <- list(1:nm, 1:nm)
    d
}

cophenetic.phylo <- function(x)
{
    n <- length(x$tip.label)
    ans <- dist.nodes(x)[1:n, 1:n]
    dimnames(ans)[1:2] <- list(x$tip.label)
    ans
}
