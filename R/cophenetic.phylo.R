## cophenetic.phylo.R (2012-08-14)

##   Pairwise Distances from a Phylogenetic Tree

## Copyright 2006-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.nodes <- function(x)
{
    x <- reorder(x) # required for the C code
    n <- Ntip(x)
    m <- x$Nnode
    nm <- n + m

    d <- .C(dist_nodes, as.integer(n), as.integer(m),
            as.integer(x$edge[, 1] - 1L), as.integer(x$edge[, 2] - 1L),
            as.double(x$edge.length), as.integer(Nedge(x)),
            double(nm * nm), NAOK = TRUE)[[7]]
    dim(d) <- c(nm, nm)
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
