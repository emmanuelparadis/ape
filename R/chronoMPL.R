## chronoMPL.R (2017-04-25)

##   Molecular Dating with Mean Path Lengths

## Copyright 2007-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

chronoMPL <- function(phy, se = TRUE, test = TRUE)
{
    if (!is.binary.phylo(phy)) stop("the tree is not dichotomous.")
    n <- length(phy$tip.label)
    m <- phy$Nnode
    N <- dim(phy$edge)[1]
    obj <- reorder(phy, "postorder")
    ndesc <- .C(node_depth, as.integer(n),
                as.integer(obj$edge[, 1]), as.integer(obj$edge[, 2]),
                as.integer(N), double(n + m), 1L)[[5]]
    s <- numeric(n + m) # sum of path lengths
    if (se) ss <- s
    if (test) Pval <- numeric(m)
    for (i in seq(1, N - 1, 2)) {
        j <- i + 1
        a <- obj$edge[i, 2]
        b <- obj$edge[j, 2]
        o <- obj$edge[i, 1]
        A <- s[a] + ndesc[a]*obj$edge.length[i]
        B <- s[b] + ndesc[b]*obj$edge.length[j]
        s[o] <- A + B
        if (se)
          ss[o] <- ss[a] + ndesc[a]^2 * obj$edge.length[i] + ss[b] +
            ndesc[b]^2 * obj$edge.length[j]
        if (test) {
            z <- abs(A/ndesc[a] - B/ndesc[b])
            tmp <- (ss[a] + ndesc[a]^2 * obj$edge.length[i])/ndesc[a]^2
            tmp <- tmp + (ss[b] + ndesc[b]^2 * obj$edge.length[j])/ndesc[b]^2
            z <- z/sqrt(tmp)
            Pval[o - n] <- 2*pnorm(z, lower.tail = FALSE)
        }
    }
    node.age <- s/ndesc
    phy$edge.length <- node.age[phy$edge[, 1]] - node.age[phy$edge[, 2]]
    if (se) attr(phy, "stderr") <- sqrt(ss[-(1:n)]/ndesc[-(1:n)]^2)
    if (test) attr(phy, "Pval") <- Pval
    phy
}
