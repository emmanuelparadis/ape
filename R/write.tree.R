## write.tree.R (2019-03-01)

##   Write Tree File in Parenthetic Format

## Copyright 2002-2019 Emmanuel Paradis, Daniel Lawson, and Klaus Schliep

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
    if (!(inherits(phy, c("phylo", "multiPhylo"))))
        stop("object \"phy\" has no trees")

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
    f.d <- paste("%.", digits, "g", sep = "")
    cp <- function(x){
        STRING[k] <<- x
        k <<- k + 1
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) cp(",")
        }
        cp(")")
        if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }

    n <- length(phy$tip.label)

    ## borrowed from phangorn:
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent))
        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

    ind <- match(1:max(phy$edge), phy$edge[, 2])

    LS <- 4*n + 5
    if (brl) LS <- LS + 4*n
    if (nodelab)  LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp(tree.prefix)
    cp("(")
    getRoot <- function(phy)
        phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) cp(",")
    }

    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}
