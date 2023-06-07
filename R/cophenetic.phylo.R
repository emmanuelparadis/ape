## cophenetic.phylo.R (2023-06-07)

##   Pairwise Distances from a Phylogenetic Tree

## Copyright 2006-2023 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.nodes <- function(x)
{
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
