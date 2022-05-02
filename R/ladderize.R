## ladderize.R (2022-04-25)

##   Ladderize a Tree

## Copyright 2007-2017 Emmanuel Paradis, 2022 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

ladderize <- function(phy, right = TRUE)
{
    desc_fun <- function(x) {
        parent <- x[, 1]
        children <- x[, 2]
        res <- vector("list", max(x))
        for (i in seq_along(parent))
            res[[parent[i]]] <- c(res[[parent[i]]], children[i])
        res
    }

    if (!is.null(phy$edge.length)) {
        el <- numeric(max(phy$edge))
        el[phy$edge[, 2]] <- phy$edge.length
    }

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    nb.edge <- dim(phy$edge)[1]

    phy <- reorder(phy, "postorder")
    N <- .C(node_depth, as.integer(nb.tip),
            as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
            as.integer(nb.edge), double(nb.tip + nb.node), 1L)[[5]]

    ii <- order(x <- phy$edge[, 1], y <- N[phy$edge[, 2]], decreasing = right)
    desc <- desc_fun(phy$edge[ii, ])

    tmp <- integer(nb.node)
    new_anc <- integer(nb.node)
    new_anc[1] <- tmp[1] <- nb.tip + 1L
    k <- nb.node
    pos <- 1L

    while (pos > 0L && k > 0) {
        current <- tmp[pos]
        new_anc[k] <- current
        k <- k - 1L
        dc <- desc[[current]]
        ind <- dc > nb.tip
        if (any(ind)) {
            l <- sum(ind)
            tmp[pos -1L + seq_len(l)] <-  dc[ind]
            pos <- pos + l - 1L
        } else {
            pos <- pos - 1L
        }
    }
    edge <- cbind(rep(new_anc, lengths(desc[new_anc])), unlist(desc[new_anc]))
    phy$edge <- edge
    if (!is.null(phy$edge.length)) phy$edge.length <- el[edge[, 2L]]
    attr(phy, "order") <- "postorder"
    phy
}
