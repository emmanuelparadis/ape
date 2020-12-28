## which.edge.R (2017-10-04)

##   Identifies Edges of a Tree

## Copyright 2004-2017 Emmanuel Paradis, 2017 Joseph W. Brown, 2017 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

getMRCA <- function(phy, tip) {
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (length(tip) < 2) return(NULL)
    Ntip <- length(phy$tip.label)
    ## <FIXME> do we need to check the value(s) in 'tip'?
    ##if (any(tip > Ntip + phy$Nnode) || any(tip < 1))
    ##    stop("value(s) out of range in 'tip'")
    ## </FIXME>
    rootnd <- Ntip + 1L

    pars <- integer(phy$Nnode) # worst case assignment, usually far too long
    tnd <- if (is.character(tip)) match(tip, phy$tip.label) else tip

    done_v <- logical(Ntip + phy$Nnode)

    ## build a lookup table to get parents faster
    pvec <- integer(Ntip + phy$Nnode)
    pvec[phy$edge[, 2]] <- phy$edge[, 1]

    ## get entire lineage for first tip
    nd <- tnd[1]
    for (k in 1:phy$Nnode) {
        nd <- pvec[nd]
        pars[k] <- nd
        if (nd == rootnd) break
    }
    pars <- pars[1:k] # delete the rest
    mrcind <- integer(max(pars))
    mrcind[pars] <- 1:k

    mrcand <- pars[1]

    ## traverse lineages for remaining tips, stop if hit common ancestor
    for (i in 2:length(tnd)) {
        cnd <- tnd[i]
        done <- done_v[cnd]
        while(!done){
            done_v[cnd] <- TRUE
            cpar <- pvec[cnd] # get immediate parent
            done <- done_v[cpar] # early exit if TRUE
            if (cpar %in% pars) {
                if (cpar == rootnd) return(rootnd) # early exit
                if(mrcind[cpar] > mrcind[mrcand]) mrcand <- cpar
                done_v[cpar] <- TRUE
                done <- TRUE
            }
            cnd <- cpar # keep going!
        }
    }
    mrcand
}

which.edge <- function(phy, group)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (is.character(group))
        group <- which(phy$tip.label %in% group)
    if (length(group) == 1)
        return(match(group, phy$edge[, 2]))

    n <- length(phy$tip.label)
    sn <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)[group]
    i <- 2L
    repeat {
        x <- unique(unlist(lapply(sn, "[", i)))
        if (length(x) != 1) break
        i <- i + 1L
    }
    d <- -(1:(i - 1L))
    x <- unique(unlist(lapply(sn, function(x) x[d])))
    match(x, phy$edge[, 2L])
}
