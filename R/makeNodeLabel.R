## makeNodeLabel.R (2009-03-22)

##   Makes Node Labels

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

# relabel derived from dist.topo.R inside .compressTipLabel,
# with an extra argument.
# makeNodeLabel is faster:
# 1. use relabel to avoid calling sort inside the loop
# 2. use digest package to avoid writing files files as with md5sum
# the node labels are the same
makeNodeLabel <- function(phy, method = "number", prefix = "Node",
                          nodeList = list(), ...)
{
    method <- sapply(method, match.arg, c("number", "md5sum", "user"),
                     USE.NAMES = FALSE)

    if ("number" %in% method)
        phy$node.label <- paste(prefix, 1:phy$Nnode, sep = "")

    if ("md5sum" %in% method) {
        phy <- relabel(phy, sort(phy$tip.label))
        nl <- character(phy$Nnode)
        pp <- prop.part(phy, check.labels = FALSE)
        labs <- attr(pp, "labels")
        for (i in seq_len(phy$Nnode)) {
            tmp <- paste0(labs[pp[[i]]], sep="\n", collapse = "")
            nl[i] <- digest(tmp, algo="md5", serialize = FALSE)
        }
        phy$node.label <- nl
    }
    if ("user" %in% method) {
        if (is.null(phy$node.label))
            phy$node.label <- character(phy$Nnode)
        nl <- names(nodeList)
        if (is.null(nl)) stop("argument 'nodeList' has no names")
        Ntip <- length(phy$tip.label)
        seq.nod <- .Call(seq_root2tip, phy$edge, Ntip, phy$Nnode)
        ## a local version to avoid the above call many times:
        .getMRCA <- function(seq.nod, tip) {
            sn <- seq.nod[tip]
            MRCA <- Ntip + 1
            i <- 2
            repeat {
                x <- unique(unlist(lapply(sn, "[", i)))
                if (length(x) != 1) break
                MRCA <- x
                i <- i + 1
            }
            MRCA
        }
        for (i in seq_along(nodeList)) {
            tips <- sapply(nodeList[[i]], grep, phy$tip.label, ...,
                           USE.NAMES = FALSE)
            j <- .getMRCA(seq.nod, unique(unlist(tips)))
            phy$node.label[j - Ntip] <- nl[i]
        }
    }
    phy
}


relabel <- function (y, ref, tip.label=TRUE) {
    label <- y$tip.label
    if (!identical(label, ref)) {
        if (length(label) != length(ref))
            stop("one tree has a different number of tips")
        ilab <- match(label, ref)
        if (any(is.na(ilab)))
            stop("one tree has different tip labels")
        ie <- match(seq_len(Ntip(y)), y$edge[, 2])
        y$edge[ie, 2] <- ilab
    }
    if(tip.label) y$tip.label <- ref
    else  y$tip.label <- NULL
    y
}



