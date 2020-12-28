## checkValidPhylo.R (2016-07-26)

##   Check the Structure of a "phylo" Object

## Copyright 2015-2016 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

checkValidPhylo <- function(phy)
{
    cat("Starting checking the validity of ", deparse(substitute(phy)), "...\n", sep = "")
    n <- m <- NULL
    if (is.null(phy$tip.label)) {
        cat("  FATAL: no element named 'tip.label' in the tree -- did you extract this tree from a \"multiPhylo\" object?\n")
    } else {
        if (!is.vector(phy$tip.label)) {
            cat("  FATAL: 'tip.label' is not a vector\n")
        } else {
            if (!is.character(phy$tip.label))
                cat("  MODERATE: 'tip.label' is not of mode \"character\"\n")
            n <- length(phy$tip.label)
            cat("Found number of tips: n =", n, "\n")
        }
    }

    if (is.null(n))
        cat("  FATAL: cannot determine the number of tips\n")

    if (is.null(phy$Nnode)) {
        cat("  FATAL: no element named 'Nnode' in the tree\n")
    } else {
        if (!is.vector(phy$Nnode)) cat("  MODERATE: 'Nnode' is not a vector\n")
        if (length(phy$Nnode) != 1) cat("  FATAL: 'Nnode' is not of length 1\n")
        if (!is.numeric(phy$Nnode)) {
            cat("  FATAL: 'Nnode' is not numeric\n")
        } else {
            if (storage.mode(phy$Nnode) != "integer")
                cat("  MODERATE: 'Nnode' is not stored as an integer\n")
        }
        if (length(phy$Nnode) == 1 && is.numeric(phy$Nnode)) {
            m <- phy$Nnode
            cat("Found number of nodes: m =", m, "\n")
        }
    }

    if (is.null(m))
        cat("  FATAL: cannot determine the number of nodes\n")

    if (is.null(phy$edge)) {
        cat("  FATAL: no element named 'edge' in the tree\n")
    } else {
        if (!is.matrix(phy$edge)) {
            cat("  FATAL: 'edge' is not a matrix\n")
        } else {
            nc <- ncol(phy$edge)
            if (nc != 2)
                cat("  FATAL: 'edge' has", nc, "columns: it MUST have 2\n")
            if (!is.numeric(phy$edge)) {
                cat("  FATAL: 'edge' is not a numeric matrix\n")
            } else {
                if (storage.mode(phy$edge) != "integer")
                    cat("  MODERATE: the matrix 'edge' is not stored as integers\n")
                if (nc == 2) {
                    if (any(phy$edge <= 0))
                        cat("  FATAL: some elements in 'edge' are negative or zero\n")
                    if (is.null(n) || is.null(m)) {
                        cat("The number of tips and/or nodes was not found: cannot check completely the 'edge' matrix\n")
                    } else {
                        tab <- tabulate(phy$edge)
                        if (length(tab) > n + m)
                            cat("  FATAL: some numbers in 'edge' are larger than 'n + m'\n")
                        if (length(tab) < n + m)
                            cat("  MODERATE: some nodes are missing in 'edge'\n")
                        if (any(tab[1:n] != 1))
                            cat("  FATAL: each tip must appear once in 'edge'\n")
                        if (any(tab[n + 1:m] < 2))
                            cat("  FATAL: all nodes should appear at least twice in 'edge'\n")
                        if (m > 1)
                            if (any(tab[n + 2:m] < 2))
                                cat("  MODERATE: some nodes are of degree 1 or less\n")
                        if (any(phy$edge[, 1] <= n & phy$edge[, 1] > 0))
                            cat("  FATAL: tips should not appear in the 1st column of 'edge'\n")
                        if (any(table(phy$edge[, 2]) > 1))
                            cat("  FATAL: nodes and tips should appear only once in the 2nd column of 'edge'\n")
                        if (any(phy$edge[, 2] == n + 1L))
                            cat("  FATAL: the root node should not appear in the 2nd column of 'edge'\n")
                    }
                }
            }
        }
    }
    cat("Done.\n")
}
