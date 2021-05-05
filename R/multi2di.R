## multi2di.R (2021-05-05)

##   Collapse or Resolve Multichotomies

## Copyright 2005-2021 Emmanuel Paradis, 2018-2021 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

multi2di <- function(phy, ...) UseMethod("multi2di")

.multi2di_ape <- function(phy, random, equiprob, n)
{
    ## n: number of tips of phy
    phy <- reorder(phy, "postorder")
    degree <- tabulate(phy$edge[, 1])
    target <- which(degree > 2)
    pos <- match(target, phy$edge[,1])
    if (!length(target)) return(phy)
    nb.edge <- dim(phy$edge)[1]
    nextnode <- n + phy$Nnode + 1L
    new.edge <- edge2delete <- NULL
    wbl <- FALSE
    if (!is.null(phy$edge.length)) {
        wbl <- TRUE
        new.edge.length <- NULL
    }

    for (i in seq_along(target)) {
        node <- target[i]
        N <- degree[node]
        ind <- pos[i] : (pos[i]+N-1L)
        desc <- phy$edge[ind, 2]
        if (random) {
          ## if we shuffle the descendants, we need to eventually
          ## reorder the corresponding branch lenghts (see below)
          ## so we store the result of sample()
            tmp <- sample(length(desc))
            desc <- desc[tmp]
            res <- if (equiprob) rtopology(N, rooted = TRUE)$edge else rtree(N)$edge
        } else {
            res <- matrix(0L, 2*N - 2, 2)
            res[, 1] <- N + rep(1:(N - 1), each = 2)
            res[, 2] <- N + rep(2:N, each = 2)
            res[seq(1, by = 2, length.out = N - 1), 2] <- 1:(N - 1)
            res[length(res)] <- N
        }
        if (wbl) {
            ## keep the branch lengths coming from `node'
            el <- numeric(dim(res)[1]) # initialized with 0's
            el[res[, 2] <= N] <-
              if (random) phy$edge.length[ind][tmp] else phy$edge.length[ind]
        }
        ## now substitute the nodes in `res'
        ## `node' stays at the "root" of these new
        ## edges whereas their "tips" are `desc'
        Nodes <- c(node, nextnode:(nextnode + N - 3L))
        res[, 1] <- Nodes[res[, 1] - N]
        tmp <- res[, 2] > N
        res[tmp, 2] <- Nodes[res[tmp, 2] - N]
        res[!tmp, 2] <- desc[res[!tmp, 2]]
        new.edge <- rbind(new.edge, res)
        edge2delete <- c(edge2delete, ind)
        if (wbl) new.edge.length <- c(new.edge.length, el)
        nextnode <- nextnode + N - 2L
        phy$Nnode <- phy$Nnode + N - 2L
    }
    phy$edge <- rbind(phy$edge[-edge2delete, ], new.edge)
    if (wbl)
        phy$edge.length <- c(phy$edge.length[-edge2delete], new.edge.length)
    if (!is.null(attr(phy, "order"))) attr(phy, "order") <- NULL
    if (!is.null(phy$node.label))
        phy$node.label <-
            c(phy$node.label, rep("", phy$Nnode - length(phy$node.label)))
    phy <- .reorder_ape(phy, "cladewise", FALSE, n, 1L) # fix by Klaus (2017-01-16)

    ## the node numbers are not in increasing order in edge[, 2]: this
    ## will confuse drop.tip and other functions (root), so renumber them
    newNb <- integer(phy$Nnode)
    newNb[1] <- n + 1L
    sndcol <- phy$edge[, 2] > n

    ## reorder node labels before changing edge:
    if (!is.null(phy$node.label)) {
        o <- 1 + rank(phy$edge[sndcol, 2])
        ## the root's label is not changed:
        phy$node.label <- phy$node.label[c(1, o)]
    }

    ## executed from right to left, so newNb is modified before phy$edge:
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2] - n] <-
        n + 2:phy$Nnode
    phy$edge[, 1] <- newNb[phy$edge[, 1] - n]
    phy
}

multi2di.phylo <- function (phy, random = TRUE, equiprob = TRUE, ...)
    .multi2di_ape(phy, random, equiprob = equiprob, length(phy$tip.label))

multi2di.multiPhylo <- function(phy, random = TRUE, equiprob = TRUE, ...)
{
    labs <- attr(phy, "TipLabel")
    oc <- oldClass(phy)
    class(phy) <- NULL
    if (is.null(labs)) phy <- lapply(phy, multi2di.phylo, random = random, equiprob = equiprob)
    else {
        phy <- lapply(phy, .multi2di_ape, random = random, equiprob = equiprob, n = length(labs))
        attr(phy, "TipLabel") <- labs
    }
    class(phy) <- oc
    phy
}

di2multi <- function(phy, ...) UseMethod("di2multi")

## by Klaus (2018-05-28)
.di2multi_ape <- function (phy, tol = 1e-08, ntips)
{
    if (is.null(phy$edge.length)) stop("the tree has no branch length")
    phy <- reorder(phy)
    e1 <- seq_len(max(phy$edge))
    ind <- which(phy$edge.length < tol & phy$edge[, 2] > ntips)
    n <- length(ind)
    if (!n) return(phy)

    for (i in ind)
        e1[phy$edge[i,2]] <-  e1[phy$edge[i,1]]

    phy$edge[, 1] <- e1[phy$edge[, 1]]
    node2del <- phy$edge[ind, 2]
    phy$edge <- phy$edge[-ind, ]
    phy$edge.length <- phy$edge.length[-ind]

    phy$Nnode <- phy$Nnode - n

    e1 <- sort(unique(phy$edge[, 1]))
    tmp <- integer(max(phy$edge))
    tmp[e1] <- ntips + seq_len(phy$Nnode)
    tmp[1:ntips] <- seq_len(ntips)

    phy$edge[] <- tmp[phy$edge]
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[-(node2del - ntips)]
    phy
}

di2multi.phylo <- function (phy, tol = 1e-08, ...)
    .di2multi_ape(phy, tol, length(phy$tip.label))

di2multi.multiPhylo <- function(phy, tol = 1e-08, ...)
{
    labs <- attr(phy, "TipLabel")
    oc <- oldClass(phy)
    class(phy) <- NULL
    if (is.null(labs)) phy <- lapply(phy, di2multi.phylo, tol = tol)
    else {
        phy <- lapply(phy, .di2multi_ape, tol = tol, ntips = length(labs))
        attr(phy, "TipLabel") <- labs
    }
    class(phy) <- oc
    phy
}
