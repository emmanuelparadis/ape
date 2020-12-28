## as.matching.R (2011-02-26)

##    Conversion Between Phylo and Matching Objects

## Copyright 2005-2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.matching <- function(x, ...) UseMethod("as.matching")

as.matching.phylo <- function(x, labels = TRUE, ...)
{
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    if (nb.tip != nb.node + 1)
        stop("the tree must be dichotomous AND rooted.")
    x <- reorder(x, "pruningwise") # cannot use "postorder" here!
    mat <- matrix(x$edge[, 2], ncol = 2, byrow = TRUE)
    nodes <- x$edge[seq(by = 2, length.out = nb.node), 1]
    ## we can use match() becoz each node appears once in `mat'
    O <- match(mat, nodes)
    new.nodes <- 1:nb.node + nb.tip
    sel <- !is.na(O)
    mat[sel] <- new.nodes[O[sel]]
    mat <- t(apply(mat, 1, sort))

    obj <- list(matching = mat)
    if (!is.null(x$edge.length))
        warning("branch lengths have been ignored")
    if (labels) {
        obj$tip.label <- x$tip.label
        if (!is.null(x$node.label))
            obj$node.label <- x$node.label[match(new.nodes, nodes)]
    }
    class(obj) <- "matching"
    obj
}

as.phylo.matching <- function(x, ...)
{
    nb.node <- dim(x$matching)[1]
    nb.tip <- nb.node + 1
    N <- 2 * nb.node
    edge <- matrix(NA, N, 2)
    new.nodes <- numeric(N + 1)
    new.nodes[N + 1] <- nb.tip + 1
    nextnode <- nb.tip + 2
    j <- 1
    for (i in nb.node:1) {
        edge[j:(j + 1), 1] <- new.nodes[i + nb.tip]
        for (k in 1:2) {
            if (x$matching[i, k] > nb.tip) {
                edge[j + k - 1, 2] <- new.nodes[x$matching[i, k]] <- nextnode
                nextnode <- nextnode + 1
            } else edge[j + k - 1, 2] <- x$matching[i, k]
        }
        j <- j + 2
    }
    obj <- list(edge = edge)
    if (!is.null(x$tip.label)) obj$tip.label <- x$tip.label
    else obj$tip.label <- as.character(1:nb.tip)
    obj$Nnode <- nb.node
    class(obj) <- "phylo"
    read.tree(text = write.tree(obj))
}
