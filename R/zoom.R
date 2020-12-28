## zoom.R (2009-07-27)

##   Zoom on a Portion of a Phylogeny

## Copyright 2003-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

zoom <- function(phy, focus, subtree = FALSE, col = rainbow, ...)
{
    if (!is.list(focus)) focus <- list(focus)
    n <- length(focus)
    for (i in 1:n)
      if (is.character(focus[[i]]))
        focus[[i]] <- which(phy$tip.label %in% focus[[i]]) # fix by Yan Wong
    if (is.function(col)) {
        col <- if (deparse(substitute(col)) == "grey") grey(1:n/n) else col(n)
    }
    ext <- vector("list", n)
    for (i in 1:n)
      ext[[i]] <- drop.tip(phy, phy$tip.label[-focus[[i]]],
                           subtree = subtree, rooted = TRUE)
    nc <- round(sqrt(n)) + 1
    nr <- ceiling(sqrt(n))
    M <- matrix(0, nr, nc)
    x <- c(rep(1, nr), 2:(n + 1))
    M[1:length(x)] <- x
    layout(M, c(1, rep(3/(nc - 1), nc - 1)))
    phy$tip.label <- rep("", length(phy$tip.label))
    colo <- rep("black", dim(phy$edge)[1])
    for (i in 1:n)
      colo[which.edge(phy, focus[[i]])] <- col[i]
    plot.phylo(phy, edge.color = colo, ...)
    for (i in 1:n)
      plot.phylo(ext[[i]], edge.color = col[i], ...)
}
