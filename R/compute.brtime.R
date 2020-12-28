## compute.brtime.R (2012-03-02)

##   Compute and Set Branching Times

## Copyright 2011-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

compute.brtime <-
    function(phy, method = "coalescent", force.positive = NULL)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    n <- length(phy$tip.label)
    m <- phy$Nnode
    N <- Nedge(phy)

    ## x: branching times (aka, node ages, depths, or heights)

    if (identical(method, "coalescent")) { # the default
        x <- 2 * rexp(m)/(as.double((m + 1):2) * as.double(m:1))
        ## x <- 2 * rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1))
        if (is.null(force.positive)) force.positive <- TRUE
    } else if (is.numeric(method)) {
        x <- as.vector(method)
        if (length(x) != m)
            stop("number of branching times given is not equal to the number of nodes")
        if (is.null(force.positive))
            force.positive <- FALSE
    }

    y <- c(rep(0, n), x) # for all nodes (terminal and internal)

    e1 <- phy$edge[, 1L] # local copies of the pointers
    e2 <- phy$edge[, 2L] #

    if (force.positive) {
        o <- .Call(seq_root2tip, phy$edge, n, m)
        list.nodes <- list(n + 1L)
        i <- 2L
        repeat {
            z <- sapply(o, "[", i)
            z <- unique(z[!(z <= n | is.na(z))])
            if (!length(z)) break
            list.nodes[[i]] <- z
            i <- i + 1L
        }
        nodes <- unlist(lapply(list.nodes, function(x) x[sample(length(x))]))
        y[nodes] <- sort(x, decreasing = TRUE)
    }

    phy$edge.length <- y[e1] - y[e2]
    phy
}
