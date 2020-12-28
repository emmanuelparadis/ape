## summary.phylo.R (2020-07-28)

##   Print Summary of a Phylogeny and "multiPhylo" operators

## Copyright 2003-2020 Emmanuel Paradis, 2006 Ben Bolker, and Klaus Schliep 2016

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

Ntip <- function(phy) UseMethod("Ntip")

Ntip.phylo <- function(phy) length(phy$tip.label)

Ntip.multiPhylo <- function(phy)
{
    labs <- attr(phy, "TipLabel")
    if (is.null(labs)) sapply(unclass(phy), Ntip.phylo)
    else setNames(rep(length(labs), length(phy)), names(phy))
}

Nnode <- function(phy, ...) UseMethod("Nnode")

Nnode.phylo <- function(phy, internal.only = TRUE, ...)
{
    if (internal.only) return(phy$Nnode)
    phy$Nnode + length(phy$tip.label)
}

Nnode.multiPhylo <- function(phy, internal.only = TRUE, ...)
{
    res <- sapply(unclass(phy), "[[", "Nnode")
    if (internal.only) return(res)
    res + Ntip.multiPhylo(phy)
}

Nedge <- function(phy) UseMethod("Nedge")

Nedge.phylo <- function(phy) dim(phy$edge)[1]

Nedge.multiPhylo <- function(phy) sapply(unclass(phy), Nedge.phylo)

summary.phylo <- function(object, ...)
{
    cat("\nPhylogenetic tree:", deparse(substitute(object)), "\n\n")
    nb.tip <- length(object$tip.label)
    nb.node <- object$Nnode
    cat("  Number of tips:", nb.tip, "\n")
    cat("  Number of nodes:", nb.node, "\n")
    if (is.null(object$edge.length))
      cat("  No branch lengths.\n")
    else {
        cat("  Branch lengths:\n")
        cat("    mean:", mean(object$edge.length), "\n")
        cat("    variance:", var(object$edge.length), "\n")
        cat("    distribution summary:\n")
        print(summary(object$edge.length)[-4])
    }
    if (is.null(object$root.edge))
      cat("  No root edge.\n")
    else
      cat("  Root edge:", object$root.edge, "\n")
    if (nb.tip <= 10) {
        cat("  Tip labels:", object$tip.label[1], "\n")
        cat(paste("             ", object$tip.label[-1]), sep = "\n")
    }
    else {
        cat("  First ten tip labels:", object$tip.label[1], "\n")
        cat(paste("                       ", object$tip.label[2:10]), sep = "\n")
    }
    if (is.null(object$node.label))
      cat("  No node labels.\n")
    else {
        if (nb.node <= 10) {
            cat("  Node labels:", object$node.label[1], "\n")
            cat(paste("              ", object$node.label[-1]), sep = "\n")
        }
        else {
            cat("  First ten node labels:", object$node.label[1], "\n")
            cat(paste("                        ", object$node.label[2:10]), sep = "\n")

        }
    }
}

### by BB:
print.phylo <- function(x, printlen = 6,...)
{
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    cat(paste("\nPhylogenetic tree with", nb.tip, "tips and", nb.node,
              "internal nodes.\n\n"))
    cat("Tip labels:\n")
    if (nb.tip > printlen) {
        cat("  ", paste(x$tip.label[1:printlen], collapse=", "), ", ...\n", sep = "")
    } else {
        cat("  ", paste(x$tip.label, collapse=", "), "\n", sep = "")
    }
    if (!is.null(x$node.label)) {
        cat("Node labels:\n")
        if (nb.node > printlen) {
            cat("  ", paste(x$node.label[1:printlen], collapse=", "), ", ...\n", sep = "")
        } else {
            cat("  ", paste(x$node.label, collapse=", "), "\n", sep = "")
        }
    }
    rlab <- if (is.rooted(x)) "Rooted" else "Unrooted"
    cat("\n", rlab, "; ", sep="")

    blen <- if (is.null(x$edge.length)) "no branch lengths." else
    "includes branch lengths."
    cat(blen, "\n", sep = "")
}

print.multiPhylo <- function(x, details = FALSE, ...)
{
    N <- length(x)
    cat(N, "phylogenetic", ifelse(N > 1, "trees\n", "tree\n"))
    if (details)
      for (i in 1:N)
        cat("tree", i, ":", length(x[[i]]$tip.label), "tips\n")
}

"[[.multiPhylo" <- function(x, i)
{
    class(x) <- NULL
    phy <- x[[i]]
    if (!is.null(attr(x, "TipLabel")))
        phy$tip.label <- attr(x, "TipLabel")
    phy
}

`$.multiPhylo` <- function(x, name) x[[name]]

"[.multiPhylo" <- function(x, i)
{
    oc <- oldClass(x)
    class(x) <- NULL
    structure(x[i], TipLabel = attr(x, "TipLabel"),
              class = oc)
}

str.multiPhylo <- function(object, ...)
{
    class(object) <- NULL
    cat('Class "multiPhylo"\n')
    str(object, ...)
}

.c_phylo_single <- function(phy) structure(list(phy), class = "multiPhylo")

c.phylo <- function(..., recursive = TRUE)
{
    obj <- list(...)
    classes <- lapply(obj, class)
    isphylo <- sapply(classes, function(x) "phylo" %in% x)
    if (all(isphylo)) {
        class(obj) <- "multiPhylo"
        return(obj)
    }
    if (!recursive) return(obj)
    ismulti <- sapply(classes, function(x) "multiPhylo" %in% x)
    if (all(isphylo | ismulti)) {
        for (i in which(isphylo)) obj[[i]] <- .c_phylo_single(obj[[i]])
        ## added by Klaus:
        for (i in which(ismulti)) obj[[i]] <- .uncompressTipLabel(obj[[i]])
        obj <- .makeMultiPhyloFromObj(obj)
    } else {
        warning('some objects not of class "phylo" or "multiPhylo": argument recursive=TRUE ignored')
    }
    obj
}

# this is an option to avoid growing the list, better check it also
# not really as important as long the list of trees is short (by Klaus)
.makeMultiPhyloFromObj <- function(obj)
{
    n <- length(obj)
    N <- lengths(obj, FALSE)
    cs <- c(0, cumsum(N))
    x <- vector("list", cs[length(cs)])
    for (i in 1:n) {
        a <- cs[i] + 1L
        b <- cs[i + 1L]
        x[a:b] <- obj[[i]]
    }
    class(x) <- "multiPhylo"
    x
}
## the original code:
##.makeMultiPhyloFromObj <- function(obj)
##{
##    n <- length(obj)
##    x <- obj[[1L]]
##    N <- length(x)
##    i <- 2L
##    while (i <= n) {
##        a <- N + 1L
##        N <- N + length(obj[[i]])
##        ## x is of class "multiPhylo", so this uses the operator below:
##        x[a:N] <- obj[[i]]
##        i <- i + 1L
##    }
##    x
##}

c.multiPhylo <- function(..., recursive = TRUE)
{
    obj <- list(...)
    if (!recursive) return(obj)
    classes <- lapply(obj, class)
    isphylo <- sapply(classes, function(x) "phylo" %in% x)
    ismulti <- sapply(classes, function(x) "multiPhylo" %in% x)
    if (!all(isphylo | ismulti)) {
        warning('some objects not of class "phylo" or "multiPhylo": argument recursive=TRUE ignored')
        return(obj)
    }
    for (i in which(isphylo)) obj[[i]] <- .c_phylo_single(obj[[i]])
    ## added by Klaus
    for (i in which(ismulti)) obj[[i]] <- .uncompressTipLabel(obj[[i]])
    .makeMultiPhyloFromObj(obj)
}

.uncompressTipLabel <- function(x)
{
    Lab <- attr(x, "TipLabel")
    if (is.null(Lab)) return(x)
    class(x) <- NULL
    for (i in 1:length(x)) x[[i]]$tip.label <- Lab
    class(x) <- "multiPhylo"
    attr(x, "TipLabel") <- NULL
    x
}

`[<-.multiPhylo` <- function(x, ..., value)
{
    ## recycling is allowed so no need to check: length(value) != length(..1)

    ## check that all elements in 'value' inherit class "phylo"
    test <- unlist(lapply(value, function(xx) !inherits(xx, "phylo")))
    if (any(test))
        stop("at least one element in 'value' is not of class \"phylo\".")

    oc <- oldClass(x)
    class(x) <- NULL

    if (is.null(attr(x, "TipLabel"))) {
        x[..1] <- value
        class(x) <- oc
        return(x)
    }

    x[..1] <- 0L # in case x needs to be elongated
    class(x) <- oc
    j <- 1L
    for (i in ..1) {
        ## x is of class "multiPhylo", so this uses the operator below:
        x[[i]] <- value[[j]]
        j <- j + 1L
    }
    x
}

`[[<-.multiPhylo` <- function(x, ..., value)
{
    if (!inherits(value, "phylo"))
        stop('trying to assign an object not of class "phylo" into an object of class "multiPhylo".')

    oc <- oldClass(x)
    class(x) <- NULL

    Lab <- attr(x, "TipLabel")

    if (!is.null(Lab)) {
        n <- length(Lab)
        if (n != length(value$tip.label))
            stop("tree with different number of tips than those in the list (which all have the same labels; maybe you want to uncompress them)")

        o <- match(value$tip.label, Lab)
        if (any(is.na(o)))
            stop("tree tip labels do not match with those in the list; maybe you want to uncompress them.")
        value$tip.label <- NULL
        ie <- match(o, value$edge[, 2])
        value$edge[ie, 2] <- 1:n
    }

    x[[..1]] <- value
    class(x) <- oc
    x
}

`$<-.multiPhylo` <- function(x, ..., value)
{
    x[[..1]] <- value
    x
}
