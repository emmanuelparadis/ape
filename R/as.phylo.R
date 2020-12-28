## as.phylo.R (2018-02-19)

##     Conversion Among Tree Objects

## Copyright 2005-2018 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

old2new.phylo <- function(phy)
{
    mode(phy$edge) <- "numeric"
    phy$Nnode <- -min(phy$edge)
    n <- length(phy$tip.label)
    NODES <- phy$edge < 0
    phy$edge[NODES] <- n - phy$edge[NODES]
    phy
}

new2old.phylo <- function(phy)
{
    NTIP <- length(phy$tip.label)
    NODES <- phy$edge > NTIP
    phy$edge[NODES] <- NTIP - phy$edge[NODES]
    mode(phy$edge) <- "character"
    phy$Nnode <- NULL
    phy
}

as.phylo <- function (x, ...)
{
    if (identical(class(x), "phylo")) return(x)
    UseMethod("as.phylo")
}

as.phylo.hclust <- function(x, ...)
{
    N <- dim(x$merge)[1]
    edge <- matrix(0L, 2*N, 2)
    edge.length <- numeric(2*N)
    ## `node' gives the number of the node for the i-th row of x$merge
    node <- integer(N)
    node[N] <- N + 2L
    cur.nod <- N + 3L
    j <- 1L
    for (i in N:1) {
        edge[j:(j + 1), 1] <- node[i]
        for (l in 1:2) {
            k <- j + l - 1L
            y <- x$merge[i, l]
            if (y > 0) {
                edge[k, 2] <- node[y] <- cur.nod
                cur.nod <- cur.nod + 1L
                edge.length[k] <- x$height[i] - x$height[y]
            } else {
                edge[k, 2] <- -y
                edge.length[k] <- x$height[i]
            }
        }
        j <- j + 2L
    }
    if (is.null(x$labels))
        x$labels <- as.character(1:(N + 1))
    obj <- list(edge = edge, edge.length = edge.length / 2,
                tip.label = x$labels, Nnode = N)
    class(obj) <- "phylo"
    reorder(obj)
}

as.phylo.phylog <- function(x, ...)
{
    tr <- read.tree(text = x$tre)
    n <- length(tr$tip.label)
    edge.length <- numeric(dim(tr$edge)[1])
    term  <- which(tr$edge[, 2] <= n)
    inte  <- which(tr$edge[, 2] > n)
    edge.length[term] <- x$leaves[tr$tip.label]
    edge.length[inte] <- x$nodes[tr$node.label][-1]
    tr$edge.length <- edge.length
    if (x$nodes["Root"] != 0) {
        tr$edge.root <- x$nodes["Root"]
        names(tr$edge.root) <- NULL
    }
    tr
}

as.hclust.phylo <- function(x, ...)
{
    if (!is.ultrametric(x)) stop("the tree is not ultrametric")
    if (!is.binary.phylo(x)) stop("the tree is not binary")
    if (!is.rooted(x)) stop("the tree is not rooted")
    n <- length(x$tip.label)
    if (n == 1) stop("needs n >= 2 observations for a classification")
    if (n == 2) {
        m <- matrix(c(-1L, -2L), 1, 2)
        bt <- x$edge.length[1]
    } else {
        x$node.label <- NULL # by Jinlong Zhang (2010-12-15)
        bt <- branching.times(x)
        N <- n - 1L
        x <- reorder(x, "postorder")
        m <- matrix(x$edge[, 2], N, 2, byrow = TRUE)
        anc <- x$edge[c(TRUE, FALSE), 1]
        bt <- bt[as.character(anc)] # 1st, reorder
        ## 2nd, sort keeping the root branching time in last (in case of
        ## rounding error if there zero-lengthed branches nead the root)
        bt <- c(sort(bt[-N]), bt[N])
        o <- match(names(bt), anc)
        m <- m[o, ]
        ## first renumber the tips:
        TIPS <- m <= n
        m[TIPS] <- -m[TIPS]
        ## then renumber the nodes:
        oldnodes <- as.numeric(names(bt))[-N]
        m[match(oldnodes, m)] <- 1:(N - 1)
        names(bt) <- NULL
    }
    obj <- list(merge = m, height = 2*bt, order = 1:n, labels = x$tip.label,
                call = match.call(), method = "unknown")
    class(obj) <- "hclust"
    obj
}

if (getRversion() >= "2.15.1") utils::globalVariables(c("network", "network.vertex.names<-"))
as.network.phylo <- function(x, directed = is.rooted(x), ...)
{
    if (is.null(x$node.label)) x <- makeNodeLabel(x)
    res <- network(x$edge, directed = directed, ...)
    network.vertex.names(res) <- c(x$tip.label, x$node.label)
    res
}

as.igraph.phylo <- function(x, directed = is.rooted(x), use.labels = TRUE, ...)
{
    ## local copy because x will be changed before evaluating is.rooted(x):
    directed <- directed
    if (use.labels) {
        if (is.null(x$node.label)) x <- makeNodeLabel(x)
        ## check added by Klaus:
        if (anyDuplicated(c(x$tip.label, x$node.label)))
            stop("Duplicated labels!")
        x$edge <- matrix(c(x$tip.label, x$node.label)[x$edge], ncol = 2)
    }
    igraph::graph_from_edgelist(x$edge, directed = directed)
}
