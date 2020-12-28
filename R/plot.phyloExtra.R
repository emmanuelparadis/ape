## plot.phyloExtra.R (2018-07-28)

##   Extra Functions for Plotting and Annotating

## Copyright 2016-2018 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plotBreakLongEdges <- function(phy, n = 1, ...) {
    o <- order(phy$edge.length, decreasing = TRUE)
    i <- o[seq_len(n)]
    phy$edge.length[i] <- max(phy$edge.length[-i])
    plot.phylo(phy, ...)
    edgelabels(edge = i, pch = 19, col = "white")
    edgelabels("//", i, frame = "n")
}

drawSupportOnEdges <- function(value, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    n <- lastPP$Ntip
    m <- lastPP$Nnode
    if (length(value) == m) value <- value[-1]
    else if (length(value) != m - 1)
        stop("incorrect number of support values")
    nodes <- 2:m + n
    i <- match(nodes, lastPP$edge[, 2])
    edgelabels(value, i, ...)
}

plotTreeTime <- function(phy, tip.dates, show.tip.label = FALSE, y.lim = NULL,
                         color = TRUE, ...)
{
    n <- Ntip(phy)
    if (length(tip.dates) != n)
        stop("number of dates does not match number of tips of the tree")
    if (is.null(y.lim)) y.lim <- c(-n/4, n)
    plot(phy, show.tip.label = show.tip.label, y.lim = y.lim, ...)
    psr <- par("usr")

    if (anyNA(tip.dates)) {
        s <- which(!is.na(tip.dates))
        tip.dates <- tip.dates[s]
    } else s <- 1:n
    range.dates <- range(as.numeric(tip.dates))

    diff.range <- range.dates[2] - range.dates[1]
    footrans <- function(x)
        psr[2] * (as.numeric(x) - range.dates[1]) / diff.range

    x1 <- footrans(tip.dates)
    y1 <- psr[3]
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    x2 <- lastPP$xx[s]
    y2 <- lastPP$yy[s]
    x1.scaled <- x1 / max(x1)
    col <-
        if (color) rgb(x1.scaled, 0, 1 - x1.scaled, alpha = .5)
        else grey(x1.scaled, alpha = 0.5)
    segments(x1, y1, x2, y2, col = col)
    at <- pretty(tip.dates)
    axis(1, at = footrans(at), labels = at)
}
