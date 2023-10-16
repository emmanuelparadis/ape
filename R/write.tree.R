## write.tree.R (2023-10-16)

##   Write Tree File in Parenthetic Format

## Copyright 2002-2023 Emmanuel Paradis, Daniel Lawson, and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

checkLabel <- function(x)
{
    ## delete all leading and trailing spaces and tabs, and
    ## the leading left and trailing right parentheses:
    ## (the syntax will work with any mix of these characters,
    ##  e.g., "    ( ( ((  " will correctly be deleted)
    x <- gsub("^[[:space:]\\(]+", "", x)
    x <- gsub("[[:space:]\\)]+$", "", x)
    ## replace all spaces and tabs by underscores:
    x <- gsub("[[:space:]]", "_", x)
    ## replace commas, colons, and semicolons with dashes:
    x <- gsub("[,:;]", "-", x)
    ## replace left and right parentheses with dashes:
    x <- gsub("[\\(\\)]", "-", x)
    x
}

write.tree <-
    function(phy, file = "", append = FALSE, digits = 10, tree.names = FALSE)
{
    if (!(inherits(phy, c("phylo", "multiPhylo"))) &&
        !all(vapply(phy, inherits, logical(1), 'phylo')))
        stop("object \"phy\" must contain trees")

    if (inherits(phy, "phylo")) phy <- c(phy)
    N <- length(phy)
    res <- character(N)

    if (is.logical(tree.names)) {
        if (tree.names) {
            tree.names <-
                if (is.null(names(phy))) character(N)
                else names(phy)
        } else tree.names <- character(N)
    }

    ## added by KS (2019-03-01):
    check_tips <- TRUE
    if (inherits(phy, "multiPhylo")) {
        if (!is.null(attr(phy, "TipLabel"))) {
            attr(phy, "TipLabel") <- checkLabel(attr(phy, "TipLabel"))
            check_tips <- FALSE
        }
    }

    ## added by EP (2019-01-23):
    phy <- .uncompressTipLabel(phy)
    class(phy) <- NULL

    for (i in 1:N)
        res[i] <- .write.tree2(phy[[i]], digits = digits,
                               tree.prefix = tree.names[i], check_tips)

    if (file == "") return(res)
    else cat(res, file = file, append = append, sep = "\n")
}

.write.tree2 <- function(phy, digits = 10, tree.prefix = "", check_tips)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
    if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste0(":%.", digits, "g")
    n <- length(phy$tip.label)

    ## terminal branches:
    terms <- match(seq_len(n), phy$edge[, 2])
    TERMS <- phy$tip.label
    if (brl) TERMS <- paste0(TERMS, sprintf(f.d, phy$edge.length[terms]))

    ## internal branches, including root edge:
    INTS <- rep(")", phy$Nnode)
    if (nodelab) INTS <- paste0(INTS, phy$node.label)
    if (brl) {
        tmp <- phy$edge.length[-terms][order(phy$edge[-terms, 2])]
        tmp <- c("", sprintf(f.d, tmp))
        if (!is.null(phy$root.edge)) tmp[1L] <- sprintf(f.d, phy$root.edge)
        INTS <- paste0(INTS, tmp)
    }

    ## find the root node:
    tmp.nodes <- unique.default(phy$edge[, 1L])
    tmp.m <- match(tmp.nodes, phy$edge[, 2L])
    root <- tmp.nodes[is.na(tmp.m)]
    if (length(root) > 1)
        stop("seems there is more than one root node")
    storage.mode(root) <- "integer"

    o <- reorderRcpp(phy$edge, n, root, 2L)
    ANC <- phy$edge[o, 1L]
    DESC <- phy$edge[o, 2L]
    NEWICK <- character(n + phy$Nnode)
    NEWICK[1:n] <- TERMS
    from <- to <- 1L
    repeat {
        thenode <- ANC[from]
        if (thenode == root) {
            to <- length(ANC)
        } else {
            while (ANC[to + 1L] == thenode) to <- to + 1L
        }
        tmp <- paste(NEWICK[DESC[from:to]], collapse = ",")
        tmp <- paste0("(", tmp, INTS[thenode - n])
        NEWICK[thenode] <- tmp
        if (thenode == root) break
        from <- to + 1L
    }
    paste0(tree.prefix, NEWICK[root], ";")
}
