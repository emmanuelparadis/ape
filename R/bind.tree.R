## bind.tree.R (2020-10-25)

##    Bind Trees

## Copyright 2003-2020 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

`+.phylo` <- function(x, y)
{
    p <- if (is.null(x$root.edge)) 0 else x$root.edge
    bind.tree(x, y, position = p)
}

bind.tree <- function(x, y, where = "root", position = 0, interactive = FALSE)
{
    nx <- length(x$tip.label)
    mx <- x$Nnode
    ROOTx <- nx + 1L
    ny <- length(y$tip.label)
    my <- y$Nnode

    if (interactive) {
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        if (lastPP$type != "phylogram" || lastPP$direction != "rightwards")
            stop("you must plot tree 'x' as a 'rightward phylogram'")
        cat("Click where you want to graft tree 'y'...\n")
        xy <- locator(1)
        d <- abs(xy$y - lastPP$yy)
        d[lastPP$xx - xy$x < 0] <- Inf
        where <- which.min(d)
        position <- lastPP$xx[where] - xy$x
        if (position < 0) position <- 0
        cat("The following parameters are used:\n")
        cat("  where =", where, " position =", position, "\n")
    } else {
        if (where == 0 || where == "root") where <- ROOTx
        if (position < 0) position <- 0
        if (where > nx + mx)
            stop("argument 'where' out of range for tree 'x'")
    }

    ## check whether both trees have branch lengths:
    switch(is.null(x$edge.length) + is.null(y$edge.length) + 1L,
           wbl <- TRUE, {
               x$edge.length <- y$edge.length <- NULL
               wbl <- FALSE
               warning("one tree has no branch lengths, they have been ignored")
           },
           wbl <- FALSE)

    yHasNoRootEdge <- is.null(y$root.edge)
    xHasNoRootEdge <- is.null(x$root.edge)

    x <- reorder(x) # fix by Veronika Boskova
    y <- reorder(y)

    ## find the row of 'where' before renumbering
    if (where == ROOTx) case <- 1 else {
        i <- which(x$edge[, 2] == where)
        case <- if (where <= nx) 2 else 3
    }
    ## case = 1 -> y is bound on the root of x
    ## case = 2 -> y is bound on a tip of x
    ## case = 3 -> y is bound on a node of x

    ## check that 'position' is correct
    if (position && wbl) {
### New in ape 3.0-1: this makes possible binding 'y' below
### a node of 'x' thus creating a new node in 'x'
###        if (!wbl)
###            stop("'position' is non-null but trees have no branch lengths")
        if (case == 1) {
            if (xHasNoRootEdge)
                stop("tree 'x' has no root edge")
            if (position > x$root.edge)
                stop("'position' is larger than x's root edge")
        } else {
            if (x$edge.length[i] < position)
                stop("'position' is larger than the branch length")
        }
    }

    ## the special case of substituting two tips:
    if (case == 2 && ny == 1 && !position) {
        x$tip.label[x$edge[i, 2]] <- y$tip.label
        if (wbl)
            x$edge.length[i] <- x$edge.length[i] + y$edge.length
        return(x)
    }

### because in all situations internal nodes need to be
### renumbered, they are changed to negatives first, and
### nodes eventually added will be numbered sequentially

    nodes <- x$edge > nx
    x$edge[nodes] <- -(x$edge[nodes] - nx) # -1, ..., -mx
    nodes <- y$edge > ny
    y$edge[nodes] <- -(y$edge[nodes] - ny + mx) # -(mx+1), ..., -(mx+my)
    ROOT <- -1L # may change later
    next.node <- -(mx + my) - 1L

    ## renumber now the tips in y:
    new.nx <- if (where <= nx && !position) nx - 1L else nx
    y$edge[!nodes] <- y$edge[!nodes] + new.nx

    ## if 'y' as a root edge, use it:
    if (!yHasNoRootEdge) {
        y$edge <- rbind(c(0, y$edge[1]), y$edge)
        ##                ^ will be filled later
        next.node <- next.node - 1L
        if (wbl) y$edge.length <- c(y$root.edge, y$edge.length)
    }

    switch(case, { # case = 1
        if (position) {
            x$root.edge <- x$root.edge - position
            x$edge <- rbind(c(next.node, x$edge[1]), x$edge)
            ROOT <- next.node
            if (wbl) x$edge.length <- c(position, x$edge.length)
        }
        if (yHasNoRootEdge) {
            j <- which(y$edge[, 1] == y$edge[1])
            y$edge[j, 1] <- ROOT
        } else y$edge[1] <- ROOT
        x$edge <- rbind(x$edge, y$edge)
        if (wbl)
            x$edge.length <- c(x$edge.length, y$edge.length)
    }, { # case = 2
        if (position) {
            x$edge[i, 2] <- next.node
            x$edge <- rbind(x$edge[1:i, ], c(next.node, where), x$edge[-(1:i), ])
            if (wbl) {
                x$edge.length[i] <- x$edge.length[i] - position
                x$edge.length <- c(x$edge.length[1:i], position, x$edge.length[-(1:i)])
            }
            i <- i + 1L
            if (yHasNoRootEdge) {
                j <- which(y$edge[, 1] == y$edge[1])
                y$edge[j, 1] <- x$edge[i, 1]
            } else y$edge[1] <- x$edge[i, 1]
        } else {
            if (yHasNoRootEdge) x$edge[i, 2] <- y$edge[1]
            else {
                ## the root edge of y is fused with the terminal edge of x
                if (wbl) y$edge.length[1] <- y$edge.length[1] + x$edge.length[i]
                y$edge[1] <- x$edge[i, 1]
                ## delete i-th edge in x:
                x$edge <- x$edge[-i, ]
                if (wbl) x$edge.length <- x$edge.length[-i]
                i <- i - 1L
            }
            x$tip.label <- x$tip.label[-where]
            ## renumber the tips that need to:
            ii <- which(x$edge[, 2] > where & x$edge[, 2] <= nx)
            x$edge[ii, 2] <- x$edge[ii, 2] - 1L
        }
        x$edge <- rbind(x$edge[1:i, ], y$edge, x$edge[-(1:i), ])
        if (wbl)
            x$edge.length <- c(x$edge.length[1:i], y$edge.length, x$edge.length[-(1:i)])
    }, { # case = 3
        if (position) {
            if (yHasNoRootEdge) {
                j <- which(y$edge[, 1] == y$edge[1])
                y$edge[j, 1] <- next.node
            } else y$edge[1] <- next.node
            x$edge <- rbind(x$edge[1:i, ], c(next.node, x$edge[i, 2]), x$edge[-(1:i), ])
            x$edge[i, 2] <- next.node
            if (wbl) {
                x$edge.length[i] <- x$edge.length[i] - position
                x$edge.length <- c(x$edge.length[1:i], position, x$edge.length[-(1:i)])
            }
            i <- i + 1L
        } else {
            if (yHasNoRootEdge) {
                j <- which(y$edge[, 1] == y$edge[1])
                y$edge[j, 1] <- x$edge[i, 2]
            } else y$edge[1] <- x$edge[i, 2]
        }
        x$edge <- rbind(x$edge[1:i, ], y$edge, x$edge[-(1:i), ])
        if (wbl)
            x$edge.length <- c(x$edge.length[1:i], y$edge.length, x$edge.length[-(1:i)])
    })

    x$tip.label <- c(x$tip.label, y$tip.label)

    if (is.null(x$node.label)) {
        if (!is.null(y$node.label))
            x$node.label <- c(rep(NA, mx), y$node.label)
    } else {
        x$node.label <-
            if (is.null(y$node.label)) c(x$node.label, rep(NA, my))
            else c(x$node.label, y$node.label)
    }

    n <- length(x$tip.label)
    x$Nnode <- dim(x$edge)[1] + 1L - n

    ## update the node labels before renumbering (this adds NA for
    ## the added nodes, and drops the label for those deleted)
    if (!is.null(x$node.label))
        x$node.label <- x$node.label[sort(-unique(x$edge[, 1]))]

    ## renumber nodes:
    newNb <- integer(x$Nnode)
    newNb[-ROOT] <- n + 1L
    sndcol <- x$edge[, 2] < 0
    ## executed from right to left, so newNb is modified before x$edge:
    x$edge[sndcol, 2] <- newNb[-x$edge[sndcol, 2]] <- n + 2:x$Nnode
    x$edge[, 1] <- newNb[-x$edge[, 1]]

    if (!is.null(x$node.label))
        x$node.label <- x$node.label[order(newNb[newNb > 0])]

    attr(x, "order") <- NULL
    reorder(x)
}
