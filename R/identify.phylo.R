## identify.phylo.R (2011-03-23)

##   Graphical Identification of Nodes and Tips

## Copyright 2008-2011 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

identify.phylo <- function(x, nodes = TRUE, tips = FALSE,
                           labels = FALSE, quiet = FALSE, ...)
{
    if (!quiet)
        cat("Click close to a node of the tree...\n")
    xy <- locator(1)
    if (is.null(xy)) return(NULL)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    ## rescale the coordinates (especially if the x- and
    ## y-scales are very different):
    pin <- par("pin")
    rescaleX <- pin[1]/max(lastPP$xx)
    xx <- rescaleX * lastPP$xx
    rescaleY <- pin[2]/max(lastPP$yy)
    yy <- rescaleY * lastPP$yy
    xy$x <- rescaleX * xy$x
    xy$y <- rescaleY * xy$y
    ## end of rescaling
    d <- (xy$x - xx)^2 + (xy$y - yy)^2 # no need to sqrt()
    NODE <- which.min(d)
    res <- list()
    if (NODE <= lastPP$Ntip) {
        res$tips <- if (labels) x$tip.label[NODE] else NODE
        return(res)
    }
    if (tips) {
        TIPS <- prop.part(x)[[NODE - lastPP$Ntip]]
        res$tips <- if (labels) x$tip.label[TIPS] else TIPS
    }
    if (nodes) {
        if (is.null(x$node.label)) labels <- FALSE
        res$nodes <- if (labels) x$node.label[NODE - lastPP$Ntip] else NODE
    }
    res
}
