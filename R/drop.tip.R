## drop.tip.R (2019-11-07)

##   Remove Tips in a Phylogenetic Tree

## Copyright 2003-2019 Emmanuel Paradis, 2017-2018 Klaus Schliep, 2018 Joseph Brown

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

keep.tip <- function(phy, tip)
{
    if (!inherits(phy, "phylo"))
        stop("object \"phy\" is not of class \"phylo\"")
    Ntip <- length(phy$tip.label)
    ## convert to indices if strings passed in
    if (is.character(tip)) {
        idx <- match(tip, phy$tip.label)
        ## stop on bad tip names
        ## alternative is to warn but proceed. not sure what stance is
        if (anyNA(idx)) {
            um <- c("umatched tip labels:\n", paste(tip[is.na(idx)], collapse = " "))
            stop(um)
        }
        tip <- idx
    } else {
        # check that passed in indices are all valid
        out.of.range <- tip > Ntip
        if (any(out.of.range)) {
            warning("some tip numbers were larger than the number of tips: they were ignored")
            tip <- tip[!out.of.range]
        }
    }
    ## get complement tip indices to drop
    toDrop <- setdiff(1:Ntip, tip)
    drop.tip(phy, toDrop)
}

extract.clade <- function(phy, node, root.edge = 0, collapse.singles = TRUE, interactive = FALSE)
{
    n <- length(phy$tip.label)
    if (interactive) {
        cat("Click close to the node...\n")
        node <- identify(phy)$nodes
    } else {
        if (length(node) > 1) {
            node <- node[1]
            warning("only the first value of 'node' has been considered")
        }
        if (is.character(node)) {
            if (is.null(phy$node.label))
                stop("the tree has no node labels")
            node <- match(node, phy$node.label) + n
            if (is.na(node)) stop("'node' not among the node labels.")
        }
        if (node <= n)
            stop("node number must be greater than the number of tips")
    }
    if (node == n + 1L) return(phy)
    keep <- prop.part(phy)[[node - n]]
    drop.tip(phy, (1:n)[-keep], root.edge = root.edge, rooted = TRUE,
             collapse.singles = collapse.singles)
}

drop.tip <-
    function(phy, tip, trim.internal = TRUE, subtree = FALSE,
             root.edge = 0, rooted = is.rooted(phy), collapse.singles = TRUE,
             interactive = FALSE)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')

    Ntip <- length(phy$tip.label)
    ## find the tips to drop:
    if (interactive) {
        cat("Left-click close to the tips you want to drop; right-click when finished...\n")
        xy <- locator()
        nToDrop <- length(xy$x)
        tip <- integer(nToDrop)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        for (i in 1:nToDrop) {
            d <- sqrt((xy$x[i] - lastPP$xx)^2 + (xy$y[i] - lastPP$yy)^2)
            tip[i] <- which.min(d)
        }
    } else {
        if (is.character(tip))
            tip <- which(phy$tip.label %in% tip)
    }

    out.of.range <- tip > Ntip
    if (any(out.of.range)) {
        warning("some tip numbers were larger than the number of tips: they were ignored")
        tip <- tip[!out.of.range]
    }

    if (!length(tip)) return(phy)

    if (length(tip) == Ntip) {
        if (Nnode(phy) < 3 || trim.internal) { # by Klaus (2018-06-21)
            warning("drop all tips of the tree: returning NULL")
            return(NULL)
        }
    }

    wbl <- !is.null(phy$edge.length)

    if (length(tip) == Ntip - 1 && trim.internal) {
        i <- which(phy$edge[, 2] == (1:Ntip)[-tip])
        res <- list(edge = matrix(2:1, 1, 2),
                    tip.label = phy$tip.label[phy$edge[i, 2]],
                    Nnode = 1L)
        class(res) <- "phylo"
        if (wbl) res$edge.length <- phy$edge.length[i]
        if (!is.null(phy$node.label))
            res$node.label <- phy$node.label[phy$edge[i, 1] - Ntip]
        return(res)
    }

    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }

    phy <- reorder(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (subtree) {
        trim.internal <- TRUE
        tr <- reorder(phy, "postorder")
        N <- .C(node_depth, as.integer(Ntip),
                as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]),
                as.integer(Nedge), double(Ntip + Nnode), 1L)[[5]]
    }

    edge1 <- phy$edge[, 1] # local copies
    edge2 <- phy$edge[, 2] #
    keep <- !logical(Nedge)

    ## delete the terminal edges given by `tip':
    keep[match(tip, edge2)] <- FALSE

    if (trim.internal) {
        ints <- edge2 > Ntip
        ## delete the internal edges that do not have anymore
        ## descendants (ie, they are in the 2nd col of `edge' but
        ## not in the 1st one)
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel)) break
            keep[sel] <- FALSE
        }
        if (subtree) {
            ## keep the subtending edge(s):
            subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
            keep[subt] <- TRUE
        }
        if (root.edge && wbl) {
            degree <- tabulate(edge1[keep])
            if (degree[ROOT] == 1) {
                j <- integer(0) # will store the indices of the edges below the new root
                repeat {
                    i <- which(edge1 == NEWROOT & keep)
                    j <- c(i, j)
                    NEWROOT <- edge2[i]
                    ## degree <- tabulate(edge1[keep]) # utile ?
                    if (degree[NEWROOT] > 1) break
                }
                keep[j] <- FALSE
                ## if (length(j) > root.edge) j <- 1:root.edge
                j <- j[1:root.edge]
                NewRootEdge <- sum(phy$edge.length[j])
                if (length(j) < root.edge && !is.null(phy$root.edge))
                    NewRootEdge <- NewRootEdge + phy$root.edge
                phy$root.edge <- NewRootEdge
            }
        }
    }

    if (!root.edge) phy$root.edge <- NULL

    ## drop the edges
    phy$edge <- phy$edge[keep, ]
    if (wbl) phy$edge.length <- phy$edge.length[keep]

    ## find the new terminal edges (works whatever 'subtree' and 'trim.internal'):
    TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])

    ## get the old No. of the nodes and tips that become tips:
    oldNo.ofNewTips <- phy$edge[TERMS, 2]

    ## in case some tips are dropped but kept because of 'subtree = TRUE':
    if (subtree) {
        i <- which(tip %in% oldNo.ofNewTips)
        if (length(i)) {
            phy$tip.label[tip[i]] <- "[1_tip]"
            tip <- tip[-i]
        }
    }

    n <- length(oldNo.ofNewTips) # the new number of tips in the tree

    ## the tips may not be sorted in increasing order in the
    ## 2nd col of edge, so no need to reorder $tip.label
    phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
    ## fix by Thomas Sibley (2017-10-28):
    if (length(tip)) phy$tip.label <- phy$tip.label[-tip]

    ## make new tip labels if necessary:
    if (subtree || !trim.internal) {
        ## get the numbers of the nodes that become tips:
        node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
        ## fix by Thomas Sibley (2017-10-28):
        new.tip.label <-
            if (!length(node2tip)) {
                character(0)
            } else if (subtree) {
                paste("[", N[node2tip], "_tips]", sep = "")
            } else {
                if (is.null(phy$node.label)) rep("NA", length(node2tip))
                else phy$node.label[node2tip - Ntip]
            }
#        if (!is.null(phy$node.label))
#            phy$node.label <- phy$node.label[-(node2tip - Ntip)]
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }

    phy$Nnode <- dim(phy$edge)[1] - n + 1L # update phy$Nnode

    ## The block below renumbers the nodes so that they conform
    ## to the "phylo" format
    newNb <- integer(Ntip + Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    storage.mode(phy$edge) <- "integer"
    if (!is.null(phy$node.label)) # update node.label if needed
        phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
    if (collapse.singles) phy <- collapse.singles(phy)
    phy
}
