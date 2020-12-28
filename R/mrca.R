## mrca.R (2017-07-28)

##   Find Most Recent Common Ancestors Between Pairs

## Copyright 2005-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

mrca <- function(phy, full = FALSE)
{
    if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo"')
    ## Get all clades:
    n <- length(phy$tip.label)
    m <- phy$Nnode
    phy <- reorder.phylo(phy, "postorder")
    BP <- bipartition2(phy$edge, n)
    N <- n + m
    ROOT <- n + 1L
    ## In the following matrix, numeric indexing will be used:
    M <- numeric(N * N)
    dim(M) <- c(N, N)

    e1 <- phy$edge[, 1]
    e2 <- phy$edge[, 2]

    ## We start at the root:
    next.node <- ROOT
    while (length(next.node)) {
        tmp <- numeric(0)
        for (anc in next.node) {
            ## Find the branches which `anc' is the ancestor...:
            id <- which(e1 == anc)
            ## ... and get their descendants:
            desc <- e2[id]
            ## `anc' is itself the MRCA of its direct descendants:
            M[anc, desc] <- M[desc, anc] <- anc
            ## Find all 2-by-2 combinations of `desc': `anc'
            ## is their MRCA:
            for (i in 1:length(desc))
                M[cbind(desc[i], desc[-i])] <- anc
            ## If one element of `desc' is a node, then the tips it
            ## leads to and the other elements of `desc' have also
            ## `anc' as MRCA!
            for (i in 1:length(desc)) {
                if (desc[i] < ROOT) next
                ## (get the tips:)
                tips <- BP[[desc[i] - n]]
                ## Same thing for the nodes...
                node.desc <- numeric(0)
                for (k in 1:m) {
                    if (k == desc[i] - n) next
                    ## If the clade of the current node is a
                    ## subset of desc[i], then it is one of its
                    ## descendants:
                    if (all(BP[[k]] %in% tips))
                      node.desc <- c(node.desc, k)
                }
                ## all nodes and tips which are descendants of
                ## `desc[i]':
                ALLDESC <- c(tips, node.desc + n)
                M[ALLDESC, desc[-i]] <- M[desc[-i], ALLDESC] <- anc
                for (j in 1:length(desc)) {
                    if (j == i || desc[j] < ROOT) next
                    tips2 <- BP[[desc[j] - n]]
                    node.desc <- numeric(0)
                    for (k in 1:m) {
                        if (k == desc[j] - n) next
                        if (all(BP[[k]] %in% tips2))
                          node.desc <- c(node.desc, k)
                    }
                    ALLDESC2 <- c(tips2, node.desc + n)
                    M[ALLDESC, ALLDESC2] <- M[ALLDESC2, ALLDESC] <- anc
                }
                ## `anc' is also the MRCA of itself and its descendants:
                M[ALLDESC, anc] <- M[anc, ALLDESC] <- anc
            }
            ## When it is done, `desc' i stored to become
            ## the new `next.node', if they are nodes:
            tmp <- c(tmp, desc[desc > n])
        }
        next.node <- tmp
    }
    M[cbind(1:N, 1:N)] <- 1:N
    if (full)
        dimnames(M)[1:2] <- list(as.character(1:N))
    else {
        M <- M[1:n, 1:n]
        dimnames(M)[1:2] <- list(phy$tip.label)
    }
    M
}
