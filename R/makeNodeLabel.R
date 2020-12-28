## makeNodeLabel.R (2009-03-22)

##   Makes Node Labels

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

makeNodeLabel <- function(phy, method = "number", prefix = "Node",
                          nodeList = list(), ...)
{
    method <- sapply(method, match.arg, c("number", "md5sum", "user"),
                     USE.NAMES = FALSE)

    if ("number" %in% method)
        phy$node.label <- paste(prefix, 1:phy$Nnode, sep = "")

    if ("md5sum" %in% method) {
        nl <- character(phy$Nnode)
        pp <- prop.part(phy, check.labels = FALSE)
        labs <- attr(pp, "labels")
        fl <- tempfile()
        for (i in seq_len(phy$Nnode)) {
            cat(sort(labs[pp[[i]]]), sep = "\n", file = fl)
            nl[i] <- tools::md5sum(fl)
        }
        unlink(fl)
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
