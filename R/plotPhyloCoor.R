## plotPhyloCoor.R (2017-05-26)

##   Coordinates of a Tree Plot

## Copyright 2008 Damien de Vienne, 2013-2017 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plotPhyloCoor <-
    function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL,
              direction = "rightwards", tip.height = NULL, ...)
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1)
        stop("found only one tip in the tree!")
    Nedge <- dim(x$edge)[1]
    if (any(tabulate(x$edge[, 1]) == 1))
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles().")
    Nnode <- x$Nnode
    if (is.null(x$edge.length)) use.edge.length <- FALSE
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    if (phyloORclado) {
        ## changed by KS:
        yy <- numeric(Ntip + Nnode)
        x <- reorder(x)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        if (!is.null(tip.height)) {
            if(!is.null(names(tip.height))) tip.height = tip.height[x$tip.label]
            yy[TIPS] <- tip.height
        } 
        else yy[TIPS] <- 1:Ntip
    }

    xe <- x$edge
    ## first reorder the tree in cladewise order to avoid cophyloplot() hanging:
    ## x <- reorder(reorder(x), order = "pruningwise") ... maybe not needed anymore (EP)
    x <- reorder(x, order = "postorder")
    ereorder <- match(x$edge[, 2], xe[, 2])

    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length)
                node.pos <- 2
        }
        if (node.pos == 1)
            yy <- .C(node_height, as.integer(x$edge[, 1]),
                     as.integer(x$edge[, 2]), as.integer(Nedge),
                     as.double(yy))[[4]]
        else {
            ans <- .C(node_height_clado, as.integer(Ntip),
                      as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                      as.integer(Nedge), double(Ntip + Nnode), as.double(yy))
            xx <- ans[[5]] - 1
            yy <- ans[[6]]
        }
        if (!use.edge.length) {
            if (node.pos != 2)
                xx <- .C(node_depth, as.integer(Ntip),
                         as.integer(x$edge[, 1]), as.integer(x$edge[, 2]),
                         as.integer(Nedge), double(Ntip + Nnode), 1L)[[5]] - 1
            xx <- max(xx) - xx
        } else {
            xx <- .C(node_depth_edgelength, as.integer(x$edge[, 1]),
                     as.integer(x$edge[, 2]), as.integer(Nedge),
                     as.double(x$edge.length), double(Ntip + Nnode))[[5]]
        }
    }

    if (phyloORclado && direction != "rightwards") {
        if (direction == "leftwards") {
            xx <- -xx
            xx <- xx - min(xx)
        }
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
            if (direction == "downwards") {
                yy <- -yy
                yy <- yy - min(yy)
            }
        }
    }
    cbind(xx, yy)
}
