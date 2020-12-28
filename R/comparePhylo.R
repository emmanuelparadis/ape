## comparePhylo.R (2018-03-13)

##   Compare Two "phylo" Objects

## Copyright 2018 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

comparePhylo <- function(x, y, plot = FALSE, force.rooted = FALSE,
                         use.edge.length = FALSE)
{
    tree1 <- deparse(substitute(x))
    tree2 <- deparse(substitute(y))
    res <- list()
    msg <- paste("=> Comparing", tree1, "with", tree2)
    res$messages <- msg
    n1 <- Ntip(x)
    n2 <- Ntip(y)
    tmp <-
        if (n1 == n2)
            paste("Both trees have the same number of tips:", n1)
        else
            paste("Trees have different numbers of tips:", n1, "and", n2)
    msg <- c(msg, tmp)
    tips1 <- x$tip.label
    tips2 <- y$tip.label
    tips12 <- match(tips1, tips2)
    tips21 <- match(tips2, tips1)
    tmp <- is.na(tips12)
    if (any(tmp))
        msg <- c(msg, paste("Tips in", tree1, "not in", tree2, ":",
                            paste(tips1[tmp], collapse = ", ")))
    tmp2 <- is.na(tips21)
    if (any(tmp2))
        msg <- c(msg, paste("Tips in", tree2, "not in", tree1, ":",
                            paste(tips2[tmp2], collapse = ", ")))
    sameTips <- FALSE
    if (!sum(tmp, tmp2)) {
        msg <- c(msg, "Both trees have the same tip labels")
        sameTips <- TRUE
    }
    m1 <- Nnode(x)
    m2 <- Nnode(y)
    tmp <-
        if (m1 == m2)
            paste("Both trees have the same number of nodes:", m1)
        else
            paste("Trees have different numbers of nodes:", m1, "and", m2)
    msg <- c(msg, tmp)
    rooted1 <- is.rooted(x)
    rooted2 <- is.rooted(y)
    tmp <- if (rooted1) {
               if (rooted2) "Both trees are rooted" else paste(tree1, "is rooted,", tree2, "is unrooted")
           } else {
               if (rooted2) paste(tree1, "is unrooted,", tree2, "is rooted") else "Both trees are unrooted"
           }
    msg <- c(msg, tmp)
    ultra1 <- ultra2 <- FALSE
    if (!is.null(x$edge.length)) ultra1 <- is.ultrametric(x)
    if (!is.null(y$edge.length)) ultra2 <- is.ultrametric(y)
    tmp <- if (ultra1) {
               if (ultra2) "Both trees are ultrametric" else paste(tree1, "is ultrametric,", tree2, "is not")
           } else {
               if (ultra2) paste(tree1, "is not ultrametric,", tree2, "is ultrametric") else "Both trees are not ultrametric"
           }
    msg <- c(msg, tmp)
    if (rooted1 && rooted2 || force.rooted) {
        key1 <- makeNodeLabel(x, "md5sum")$node.label
        key2 <- makeNodeLabel(y, "md5sum")$node.label
        mk12 <- match(key1, key2)
        mk21 <- match(key2, key1)
        if (any(tmp <- is.na(mk12))) {
            nk <- sum(tmp)
            msg <- c(msg, paste(nk, if (nk == 1) "clade" else "clades", "in", tree1, "not in", tree2))
        }
        if (plot) {
            layout(matrix(1:2, 1, 2))
            plot(x, use.edge.length = use.edge.length, main = tree1)
            nodelabels(node = which(tmp) + n1, pch = 19, col = "blue", cex = 2)
            legend("topleft", legend = paste("Clade absent in", tree2), pch = 19, col = "blue")
        }
        if (any(tmp <- is.na(mk21))) {
            nk <- sum(tmp)
            msg <- c(msg, paste(nk, if (nk == 1) "clade" else "clades", "in", tree2, "not in", tree1))
        }
        if (plot) {
            plot(y, use.edge.length = use.edge.length, main = tree2)
            nodelabels(node = which(tmp) + n2, pch = 19, col = "red", cex = 2)
            legend("topleft", legend = paste("Clade absent in", tree1), pch = 19, col = "red")
        }
        nodes1 <- which(!is.na(mk12))
        nodes2 <- mk12[!is.na(mk12)]
        if (ultra1 && ultra2) {
            bt1 <- branching.times(x)
            bt2 <- branching.times(y)
            BT <- data.frame(paste0(bt1[nodes1], " (", nodes1 + n1, ")"),
                             paste0(bt2[nodes2], " (", nodes2 + n2, ")"))
            names(BT) <- c(tree1, tree2)
            res$BT <- BT
            msg <- c(msg, "Branching times of clades in common between both trees: see ..$BT
(node number in parentheses)")
        }
        if (!is.null(nl1 <- x$node.label) && !is.null(nl2 <- y$node.label)) {
            NODES <- data.frame(paste0(nl1[nodes1], " (", nodes1 + n1, ")"),
                                paste0(nl2[nodes2], " (", nodes2 + n2, ")"))
            names(NODES) <- c(tree1, tree2)
            res$NODES <- NODES
            msg <- c(msg, "Node labels of clades in common between both trees: see ..$NODES
(node number in parentheses)")
        }
    }
    if (!force.rooted && !rooted1 && !rooted2 && sameTips && m1 == m2) {
        TR <- .compressTipLabel(c(x, y))
        bs <- bitsplits(TR)
        common.splits <- which(bs$freq == 2L)
        ncs <- length(common.splits)
        tmp <-
            if (ncs)
                paste(ncs, if (ncs == 1) "split" else "splits", "in common")
            else "No split in common"
        msg <- c(msg, tmp)
        if (plot) {
            co <- "black"#rgb(0, 0, 1, 0.7)
            layout(matrix(1:2, 1, 2))
            edgecol1 <- rep("black", Nedge(x))
            edgew1 <- rep(1, Nedge(x))
            edgecol2 <- rep("black", Nedge(y))
            edgew2 <- rep(1, Nedge(y))
            if (ncs) {
                f <- function(x)
                    unlist(lapply(ONEwise(x), paste, collapse = "\r"))
                ##pp <-  f(as.prop.part(bs, include.trivial = TRUE))
                pp <-  as.prop.part(bs)
                pp1 <- f(prop.part(TR[[1]]))
                pp2 <- f(prop.part(TR[[2]]))
                one2n <- 1:n1
                for (i in common.splits) {
                    p <- pp[[i]]
                    split <- paste(p, collapse = "\r")
                    k1 <- match(split, pp1)
                    k2 <- match(split, pp2)
                    if (!length(k1)) {
                        split <- paste(one2n[-p], collapse = "\r")
                        k1 <- match(split, pp1)
                        if (!length(k2)) k2 <- match(split, pp2)
                    }
                    e1 <- match(k1 + n1, TR[[1]]$edge[, 2])
                    e2 <- match(k2 + n2, TR[[2]]$edge[, 2])
                    edgecol1[e1] <- edgecol2[e2] <- co
                    edgew1[e1] <- edgew2[e2] <- 5
                }
            }
            plot(TR[[1]], "u", use.edge.length = use.edge.length,
                 edge.color = edgecol1, edge.width = edgew1, main = tree1, cex = 1.3, font =1)
            legend("bottomright", legend = "Split present in both trees",
                   lty = 1, col = "black", lwd = 5)
            plot(TR[[2]], "u", use.edge.length = use.edge.length,
                 edge.color = edgecol2, edge.width = edgew2, main = tree2, cex = 1.3, font =1)
        }
    }
    res$messages <- paste0(msg, ".")
    class(res) <- "comparePhylo"
    res
}

print.comparePhylo <- function(x, ...)
{
    cat(x$messages, sep = "\n")
    cat("\n")
    x$messages <- class(x) <- NULL
    if (length(x)) print.default(x)
}
