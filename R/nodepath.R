## nodepath.R (2014-11-06)

##   Find Paths of Nodes

## Copyright 2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

nodepath <- function(phy, from = NULL, to = NULL)
{
    if (!inherits(phy, "phylo"))
        stop("object \"phy\" is not of class \"phylo\"")
    n <- length(phy$tip.label)
    m <- phy$Nnode
    root2tip <- .Call(seq_root2tip, phy$edge, n, m)
    if (is.null(from) || is.null(to)) return(root2tip)
    if (from < 1 || from > n + m) stop("'from' out of range")
    if (to < 1 || to > n + m) stop("'to' out of range")
    if (from == to) return(to)

    ## find the first occurrence of 'x' in the list root2tip
    foo <- function(x) {
        if (x <= n) return(x) # if x is a tip
        if (x == n + 1L) return(1L) # if x is the root
        i <- 1L
        repeat {
            if (any(root2tip[[i]] == x)) break
            i <- i + 1L
        }
        i
    }

    i <- foo(from)
    j <- foo(to)

    ## find path of nodes in a single vector 'seq' from root2tip
    findPath <- function(from, to, seq) {
        i <- which(seq == from)
        j <- which(seq == to)
        seq[i:j]
    }

    if (i == j) return(findPath(from, to, root2tip[[i]]))

    ## find the MRCA of 'from' and 'to'
    A <- root2tip[[i]]
    B <- root2tip[[j]]
    MRCA <- n + 1L # start from the root
    k <- 2L
    repeat {
        if (A[k] != B[k]) break
        MRCA <- A[k]
        k <- k + 1L
    }

    x <- findPath(MRCA, from, A)
    y <- findPath(MRCA, to, B)
    c(rev(x), y[-1])
}
