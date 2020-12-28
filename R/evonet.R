## evonet.R (2017-07-28)

##   Evolutionary Networks

## Copyright 2011-2012 Emmanuel Paradis, 2017 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

evonet <- function(phy, from, to = NULL)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo".')
    if (!is.rooted(phy))
        warning("the tree is unrooted")
    x <- phy

    if (is.null(to)) {
        if (is.data.frame(from))
            from <- as.matrix(from)
        if (!is.matrix(from))
            stop("'from' must be a matrix or a data frame if 'to' is not given")
        if (ncol(from) > 2) {
            warning("'from' has more than two columns: only the first two will be used.")
            ret <- from[, 1:2]
        } else if (ncol(from) < 2) {
            stop("'from' must have at least two columns")
        } else ret <- from
    } else {
        from <- as.vector(from)
        to <- as.vector(to)
        if (length(from) != length(to))
            stop("'from' and 'to' not of the same length after coercing as vectors")
        ret <- cbind(from, to)
    }

    ## check that values are not out of range:
    storage.mode(ret) <- "integer"
    if (any(is.na(ret)))
        stop("some values are NA's after coercing as integers")
    if (any(ret < 0) || any(ret > Ntip(phy) + phy$Nnode))
        stop("some values are out of range")

    x$reticulation <- ret
    class(x) <- c("evonet", "phylo")
    x
}

as.phylo.evonet <- function(x, ...)
{
    x$reticulation <- NULL
    class(x) <- "phylo"
    x
}

plot.evonet <- function(x, col = "blue", lty = 1, lwd = 1, alpha = 0.5,
                        arrows = 0, arrow.type = "classical", ...)
{
    ## changed 5/24/17 by Klaus
    plot.phylo(x, ...)
    edges(x$reticulation[, 1], x$reticulation[, 2],
          col = rgb(t(col2rgb(col)), alpha = 255 * alpha,
                    maxColorValue = 255),
          lty = lty, lwd = lwd, arrows = arrows, type = arrow.type)
}

as.networx.evonet <- function(x, weight = NA, ...)
{
    if (any(x$reticulation <= Ntip(x)))
        stop("some tips are involved in reticulations: cannot convert to \"networx\"")
    x <- reorder(x, "postorder")
    ned <- Nedge(x)
    nrt <- nrow(x$reticulation)
    x$edge <- rbind(x$edge, x$reticulation)
    colnames(x$edge) <- c("oldNodes", "newNodes")
    x$reticulation <- NULL
    x$edge.length <- c(x$edge.length, rep(weight, length.out = nrt))
    x$split <- c(1:ned, 1:nrt)
    class(x) <- c("networx", "phylo")
    x
}

as.network.evonet <- function(x, directed = TRUE, ...)
{
    class(x) <- NULL
    x$edge <- rbind(x$edge, x$reticulation)
    as.network.phylo(x, directed = directed, ...)
}

as.igraph.evonet <- function(x, directed = TRUE, use.labels = TRUE, ...)
{
    class(x) <- NULL
    x$edge <- rbind(x$edge, x$reticulation)
    ## added check by Klaus (2017-05-26)
    if (use.labels) {
        if (!is.null(x$node.label)){
            tmp <- nchar(x$node.label)
            if (any(tmp == 0)){
                newLabel <- paste0("number", 1:x$Nnode)
                x$node.label[tmp == 0] <- newLabel[tmp == 0]
            }
        }
        if (any(duplicated(c(x$tip.label, x$node.label))))
            stop("Duplicated labels!")
    }
    as.igraph.phylo(x, directed = directed, use.labels = use.labels, ...)
}

print.evonet <- function(x, ...)
{
    nr <- nrow(x$reticulation)
    cat("\n    Evolutionary network with", nr, "reticulation")
    if (nr > 1) cat("s")
    cat("\n\n               --- Base tree ---")
    print.phylo(as.phylo(x))
}

## new stuff by Klaus (2017-05-26)

reorder.evonet <- function(x, order = "cladewise", index.only = FALSE, ...)
{
    reticulation <- x$reticulation
    y <- reorder(as.phylo(x), order = order, index.only = index.only, ...)
    if (index.only) return(y)
    y$reticulation <- reticulation
    class(y) <- c("evonet", "phylo")
    y
}

## requires topo_sort from igraph, behaviour different from phylo
## (postorder seems to work fine)
## if no singletons are in edge reorder.phylo could be used
## if (getRversion() >= "2.15.1") utils::globalVariables(c("topo_sort", "graph"))
## reorder.evonet <- function(x, order = "cladewise", index.only = FALSE, ...)
## {
##     order <- match.arg(order, c("cladewise", "postorder"))
##     if (!is.null(attr(x, "order")))
##         if (attr(x, "order") == order) return(x)
##     g <- graph(t(x$edge))
##     neword <- if (order == "cladewise") topo_sort(g, "out") else topo_sort(g, "in")
##     neworder <- order(match(x$edge[, 1], neword))
##     if (index.only) return(neworder)
##     x$edge <- x$edge[neworder, ]
##     if (!is.null(x$edge.length)) x$edge.length <- x$edge.length[neworder]
##     attr(x, "order") <- order
##     x
## }

as.evonet <- function(x, ...)
{
    if (inherits(x, "evonet")) return(x)
    UseMethod("as.evonet")
}

as.evonet.phylo <- function(x, ...)
{
    pos <- grep("#", x$tip.label)
    ind <- match(pos, x$edge[, 2])
    reticulation <- x$edge[ind, , drop = FALSE]
    edge <- x$edge[-ind, , drop = FALSE]
    nTips <- as.integer(length(x$tip.label))
    reticulation[, 2] <- as.integer(match(x$tip.label[pos], x$node.label) + nTips)
    for (i in sort(pos, TRUE)) {
        edge[edge > i ] <- edge[edge > i] - 1L
        reticulation[reticulation > i] <- reticulation[reticulation > i] - 1L
    }
    x$edge <- edge
    x$reticulation <- reticulation
    if (!is.null(x$edge.length)) x$edge.length <- x$edge.length[-ind]
    x$tip.label <- x$tip.label[-pos]
    class(x) <- c("evonet", "phylo")
    x
}

## requires new version of clado.build and tree.build
read.evonet <- function(file = "", text = NULL, comment.char = "", ...)
{
    x <- read.tree(file = file, text = text, comment.char = comment.char, ...)
    as.evonet.phylo(x)
}

.evonet2phylo <-  function(x)
{
    nTips <- as.integer(length(x$tip.label))
    if (!is.null(x$edge.length)) {
        nd <- node.depth.edgelength(x)
        x$edge.length <- c(x$edge.length, nd[x$reticulation[, 2]] - nd[x$reticulation[, 1]])
    }
    if (!is.null(x$node.label))
        x$tip.label <- c(x$tip.label, x$node.label[x$reticulation[, 2] - nTips])
    else {
        newLabels <- paste0("#H", x$reticulation[, 2])
        x$tip.label <- c(x$tip.label, newLabels)
        x$node.label <- rep("", x$Nnode)
        ind <- which((x$reticulation[, 2] > nTips) & !duplicated(x$reticulation[, 2]))
        x$node.label[x$reticulation[ind, 2] - nTips] <- newLabels[ind]
    }
    nrets <- as.integer(nrow(x$reticulation))
    x$edge[x$edge > nTips] <-  x$edge[x$edge > nTips] + nrets
    x$reticulation[, 1] <- x$reticulation[, 1] + nrets
    x$reticulation[, 2] <- nTips + (1L:nrets)
    x$edge <- rbind(x$edge, x$reticulation)
    x$reticulation <- NULL
    attr(x, "order") <- NULL
    class(x) <- "phylo"
    x
}

write.evonet <- function(x, file = "", ...)
{
    x <- .evonet2phylo(x)
    write.tree(x, file = file, ...)
}

Nedge.evonet <- function(phy) dim(phy$edge)[1] + dim(phy$reticulation)[1]
