## rtree.R (2020-11-14)

##   Generates Trees

## Copyright 2004-2020 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.N <- unname(howmanytrees(792, labeled = FALSE, detail = TRUE))
.lN <- log10(.N)
.xi <- 2.477993
.lxi <- log10(.xi)
.log10sumxipow0to9 <- log10(sum(.xi^(0:9)))
.getProb4rtree <- function(n) {
    a <- 1:floor(n/2)
    b <- n - a
    x <- .N[a] * .N[b]
    x <- x / max(x)
    p <- x / sum(x)
    if (all(is.finite(p))) return(p)
    ## use the log10-scale
    foo <- function(n) 0.3941 * n - 4.153
    lNa <- .lN[a]
    lNa[a > 792] <- foo(a[a > 792])
    lNb <- .lN[b]
    lNb[b > 792] <- foo(b[b > 792])
    lx <- lNa + lNb # log10-scale
    ## we use the first 10 terms of the approximation (see the vignette)
    sx <- (n - 1) * .lxi + (n - 10) * .lxi + .log10sumxipow0to9
    p <- lx - sx
    p <- p - min(p) # rescale
    p / sum(p)
}
rtree <- function(n, rooted = TRUE, tip.label = NULL, br = runif,
                  equiprob = FALSE, ...)
{
    ## as.integer(runif()) is more efficient than sample.int() but we
    ## have to keep the latter for the default of 'equiprob' because
    ## this is used in other packages with set.seed()
    if (equiprob) {
        bar <- function(n) {
            if (n < 4L) return(1L)
            p <- .getProb4rtree(n)
            sample.int(floor(n / 2), 1L, FALSE, p, FALSE)
            ##if (n < 4L) return(1L)
            ##as.integer(runif(1L, 0, floor(n / 2))) + 1L
        }
    } else {
        bar <- function(n) sample.int(n - 1L, 1L, FALSE, NULL, FALSE)
    }

    foo <- function(n, pos) {
        n1 <- bar(n)
        n2 <- n - n1
        po2 <- pos + 2L * n1 - 1L
        edge[c(pos, po2), 1L] <<- nod
        nod <<- nod + 1L
        if (n1 > 2L) {
            edge[pos, 2L] <<- nod
            foo(n1, pos + 1L)
        } else if (n1 == 2L) {
            edge[pos + 1:2, 1L] <<- edge[pos, 2L] <<- nod
            nod <<- nod + 1L
        }
        if (n2 > 2L) {
            edge[po2, 2L] <<- nod
            foo(n2, po2 + 1L)
        } else if (n2 == 2L) {
            edge[po2 + 1:2, 1L] <<- edge[po2, 2L] <<- nod
            nod <<- nod + 1L
        }
    }

    if (n < 1) stop("a tree must have at least 1 tip")
    if (n < 3 && !rooted) stop("an unrooted tree must have at least 3 tips")
    n <- as.integer(n)

    ## make the tip labels:
    if (is.null(tip.label)) {
        tip.label <- paste0("t", 1:n)
    } else {
        tip.label <- as.character(tip.label)
        Nlabs <- length(tip.label)
        if (!Nlabs) {
            warning("vector 'tip.label' of length zero: generating tip labels")
            tip.label <- paste0("t", seq_len(n))
        } else if (Nlabs > n) {
            warning("vector 'tip.label' longer than 'n': was shorten")
            tip.label <- tip.label[1:n]
        } else if (Nlabs < n) {
            warning("vector 'tip.label' shorter than 'n': was recycled")
            tip.label <- rep(tip.label, length.out = n)
        }
    }

    if (n == 1L) { # rooted case with n = 1
        nbr <- 1L
        edge <- matrix(2:1, 1L, 2L)
    } else { # all other cases
        nbr <- 2L * n - 3L + rooted
        edge <- matrix(NA_integer_, nbr, 2L)
    }

    if (rooted) {
        if (n == 2L) {
            edge[] <- c(3L, 3L, 1L, 2L)
        } else if (n == 3L) {
            edge[] <- c(4L, 5L, 5L, 4L, 5L, 1:3)
        } else if (n > 3L) {
            nod <- n + 1L
            foo(n, 1L)
            ## slightly more efficient than affecting the tip numbers in foo():
            i <- which(is.na(edge[, 2L]))
            edge[i, 2L] <- 1:n
        }
    } else { # unrooted case
        if (n == 3L) {
            edge[] <- c(4L, 4L, 4L, 1:3)
        } else if (n == 4L) {
            edge[] <- c(5L, 6L, 6L, 5L, 5L, 6L, 1:4)
        } else if (n == 5L) {
            edge[] <- c(6L, 6L, 6L, 7L, 7L, 8L, 8L, 1L, 2L, 7L, 3L, 8L, 4L, 5L)
        } else { # n > 5
            ## generate a rooted tree without branch lengths and unroot it
            phy <- rtree(n, tip.label = tip.label, br = NULL, equiprob = equiprob, ...)
            phy <- .unroot_ape(phy, n)
        }
    }
    if (!exists("phy", inherits = FALSE)) {
        phy <- list(edge = edge, tip.label = sample(tip.label))
        phy$Nnode <- if (n == 1L) 1L else n - 2L + as.integer(rooted)
        class(phy) <- "phylo"
        attr(phy, "order") <- "cladewise"
    }
    if (!is.null(br)) {
        phy$edge.length <-
            if (is.function(br)) br(nbr, ...) else rep(br, length.out = nbr)
    }
    phy
}

rcoal <- function(n, tip.label = NULL, br = "coalescent", ...)
{
    n <- as.integer(n)
    nbr <- 2*n - 2
    edge <- matrix(NA, nbr, 2)
    ## coalescence times by default:
    x <- if (is.character(br)) 2*rexp(n - 1)/(as.double(n:2) * as.double((n - 1):1))
    else if (is.numeric(br)) rep(br, length.out = n - 1) else br(n - 1, ...)
    if (n == 2) {
        edge[] <- c(3L, 3L, 1:2)
        edge.length <- rep(x, 2)
    } else if (n == 3) {
        edge[] <- c(4L, 5L, 5L, 4L, 5L, 1:3)
        edge.length <- c(x[c(2, 1, 1)], sum(x))
    } else {
        edge.length <- numeric(nbr)
        h <- numeric(2*n - 1)
        node.height <- cumsum(x)
        pool <- 1:n
        nextnode <- 2L*n - 1L
        for (i in 1:(n - 1)) {
            y <- sample(pool, size = 2)
            ind <- (i - 1)*2 + 1:2
            edge[ind, 2] <- y
            edge[ind, 1] <- nextnode
            edge.length[ind] <- node.height[i] - h[y]
            h[nextnode] <- node.height[i]
            pool <- c(pool[! pool %in% y], nextnode)
            nextnode <- nextnode - 1L
        }
    }
    phy <- list(edge = edge, edge.length = edge.length)
    if (is.null(tip.label))
        tip.label <- paste("t", 1:n, sep = "")
    phy$tip.label <- sample(tip.label)
    phy$Nnode <- n - 1L
    class(phy) <- "phylo"
    phy <- reorder(phy)
    ## to avoid crossings when converting with as.hclust:
    phy$edge[phy$edge[, 2] <= n, 2] <- 1:n
    phy
}

rmtree <- function(N, n, rooted = TRUE, tip.label = NULL, br = runif, equiprob = FALSE, ...)
{
    a <- replicate(N, rtree(n, rooted = rooted, tip.label =  tip.label,
                            br = br, equiprob = equiprob, ...),
                   simplify = FALSE)
    class(a) <- "multiPhylo"
    a
}

stree <- function(n, type = "star", tip.label = NULL)
{
    type <- match.arg(type, c("star", "balanced", "left", "right"))
    n <- as.integer(n)
    if (type == "star") {
        N <- n
        m <- 1L
    } else {
        m <- n - 1L
        N <- n + m - 1L
    }
    edge <- matrix(0L, N, 2)

    switch(type, "star" = {
        edge[, 1] <- n + 1L
        edge[, 2] <- 1:n
    }, "balanced" = {
        if (log2(n) %% 1)
            stop("'n' is not a power of 2: cannot make a balanced tree")
        foo <- function(node, size) {
            if (size == 2) {
                edge[c(i, i + 1L), 1L] <<- node
                edge[c(i, i + 1L), 2L] <<- c(nexttip, nexttip + 1L)
                nexttip <<- nexttip + 2L
                i <<- i + 2L
            } else {
                for (k in 1:2) { # do the 2 subclades
                    edge[i, ] <<- c(node, nextnode)
                    nextnode <<- nextnode + 1L
                    i <<- i + 1L
                    foo(nextnode - 1L, size/2)
                }
            }
        }
        i <- 1L
        nexttip <- 1L
        nextnode <- n + 2L
        foo(n + 1L, n)
    }, "left" = {
        edge[c(seq.int(from = 1, to = N - 1, by = 2), N), 2L] <- 1:n
        nodes <- (n + 1L):(n + m)
        edge[seq.int(from = 2, to = N - 1, by = 2), 2L] <- nodes[-1]
        edge[, 1L] <- rep(nodes, each = 2)
    }, "right" = {
        nodes <- (n + 1L):(n + m)
        edge[, 1L] <- c(nodes, rev(nodes))
        edge[, 2L] <- c(nodes[-1], 1:n)
    })

    if (is.null(tip.label))
        tip.label <- paste("t", 1:n, sep = "")
    phy <- list(edge = edge, tip.label = tip.label, Nnode = m)
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

.check.tip.label <- function(tip.label, n, prefix = "t")
{
    if (is.null(tip.label)) return(paste0(prefix, seq_len(n)))
    tip.label <- as.character(tip.label)
    Nlabs <- length(tip.label)
    if (!Nlabs) {
        warning("vector 'tip.label' of length zero: generating tip labels",
                call. = FALSE)
        return(paste0(prefix, seq_len(n)))
    }
    if (Nlabs > n) {
        warning("vector 'tip.label' longer than 'n': was shorten",
                call. = FALSE)
        return(tip.label[1:n])
    }
    if (Nlabs < n) {
        warning("vector 'tip.label' shorter than 'n': was recycled",
                call. = FALSE)
        return(rep(tip.label, length.out = n))
    }
    tip.label
}

rtopology <- function(n, rooted = FALSE, tip.label = NULL, br = runif, ...)
{
    n <- as.integer(n)
    if (n < 1)
        stop("a tree must have at least 1 tip")
    if (n < 3 && !rooted)
        stop("an unrooted tree must have at least 3 tips")
    if (n < 4)
        return(rtree(n, rooted = rooted, tip.label = tip.label, br = br, ...))
    nb <- n - 3L
    x <- as.integer(runif(nb) * seq(3, by = 2, length.out = nb)) + 1L

    tip.label <- .check.tip.label(tip.label, n)
    Nnode <- n - 2L

    TIPS <- sample.int(n) # permute the tips beforehand
    N <- 3L * n - 6L
    edge <- matrix(NA_integer_, N, 2L)
    alive <- logical(N)
    alive[1:3] <- TRUE
    Nalive <- 3L
    e <- 1:3
    ROOT <- n + 1L
    edge[1:3] <- ROOT
    nextnode <- ROOT + 1L
    edge[1:3 + N] <- TIPS[1:3]
    i <- 4L

    while (i <= n) {
        ## draw a branch among the alive ones
        k <- which(alive)[x[i - 3L]] # find its location in edge
        alive[k] <- FALSE # delete 1 branch
        e <- e + 3L # add 3 new branches
        alive[e] <- TRUE
        edge[e[1]] <- edge[k]
        edge[e[1] + N] <- nextnode
        edge[e[2:3]] <- nextnode
        edge[e[2] + N] <- edge[k + N]
        edge[e[3] + N] <- TIPS[i]
        nextnode <- nextnode + 1L
        Nalive <- Nalive + 2L
        i <- i + 1L
    }

    edge <- edge[alive, ]
    phy <- list(edge = edge, tip.label = tip.label, Nnode = Nnode)
    class(phy) <- "phylo"
    phy <- reorder(phy)

    if (rooted) {
        ## exclude the root partition and add the terminal trivial
        ## partitions
        og <- sample.int(n + Nnode - 1L, 1L)
        if (og > n) {
            pp <- prop.part(phy)[-1L]
            og <- pp[[og - n]]
        }
        phy <- root.phylo(phy, og, resolve.root = TRUE)
    }
    nbr <- Nedge.phylo(phy)
    if (!is.null(br)) {
        phy$edge.length <- if (is.function(br))
            br(nbr, ...)
        else rep(br, length.out = nbr)
    }
    phy
}

rmtopology <- function(N, n, rooted = FALSE, tip.label = NULL, br = runif, ...)
{
    a <- replicate(N, rtopology(n, rooted = rooted, tip.label = tip.label, br = br, ...), simplify = FALSE)
    class(a) <- "multiPhylo"
    a
}
