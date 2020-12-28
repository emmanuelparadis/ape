## MPR.R (2010-08-10)

##   Most Parsimonious Reconstruction

## Copyright 2010 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

MPR <- function(x, phy, outgroup)
{
    if (is.rooted(phy))
        stop("the tree must be unrooted")
    if (!is.binary.phylo(phy))
        stop("the tree must be fully dichotomous")
    if (length(outgroup) > 1L)
        stop("outgroup must be a single tip")
    if (is.character(outgroup))
        outgroup <- which(phy$tip.label == outgroup)

    if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label))
            x <- x[phy$tip.label]
        else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
    }

    n <- length(phy$tip.label)
    if (is.null(phy$node.label))
        phy$node.label <- n + 1:(phy$Nnode)

    phy <- drop.tip(root(phy, outgroup), outgroup)
    n <- n - 1L
    m <- phy$Nnode
    phy <- reorder(phy, "postorder")

    root.state <- x[outgroup]
    I <- as.integer(x[-outgroup])
    I[n + 1:m] <- NA
    I <- cbind(I, I) # interval map

    med <- function(x) {
        i <- length(x)/2
        sort(x)[c(i, i + 1L)]
    }

    ## 1st pass
    s <- seq(from = 1, by = 2, length.out = m)
    anc <- phy$edge[s, 1]
    des <- matrix(phy$edge[, 2], ncol = 2, byrow = TRUE)
    for (i in 1:m) I[anc[i], ] <- med(I[des[i, ], ])

    ## 2nd pass
    out <- matrix(NA, m, 2)
    colnames(out) <- c("lower", "upper")
    ## do the most basal node before looping
    Iw <- as.vector(I[des[m, ], ]) # interval maps of the descendants
    out[anc[m] - n, ] <- range(med(c(root.state, root.state, Iw)))
    for (i in (m - 1):1) {
        j <- anc[i]
        Iw <- as.vector(I[des[i, ], ]) # interval maps of the descendants
        k <- which(phy$edge[, 2] == j) # find the ancestor
        tmp <- out[phy$edge[k, 1] - n, ]
        out[j - n, 1] <- min(med(c(tmp[1], tmp[1], Iw)))
        out[j - n, 2] <- max(med(c(tmp[2], tmp[2], Iw)))
    }
    rownames(out) <- phy$node.label
    out
}
