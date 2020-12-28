## dist.topo.R (2020-07-20)

##      Topological Distances, Tree Bipartitions,
##   Consensus Trees, and Bootstrapping Phylogenies

## Copyright 2005-2020 Emmanuel Paradis, 2016-2017 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.getTreesFromDotdotdot <- function(...)
{
    obj <- list(...)
    if (length(obj) == 1 && !inherits(obj[[1]], "phylo")) obj <- obj[[1]]
    obj
}

dist.topo <- function(x, y = NULL, method = "PH85")
{
    method <- match.arg(method, c("PH85", "score"))
    if (!is.null(y)) x <- c(x, y)
    testroot <- any(is.rooted(x))
    n <- length(x) # number of trees
    res <- numeric(n * (n - 1) /2)
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("tree", 1:n)

    if (method == "PH85") {
        if (testroot)
            warning("Some trees were rooted: topological distances may be spurious.")

        x <- .compressTipLabel(x)
        ntip <- length(attr(x, "TipLabel"))
        nnode <- sapply(x, Nnode)

        foo <- function(phy, ntip) {
            phy <- reorder(phy, "postorder")
            ans <- ONEwise(bipartition2(phy$edge, ntip))
            sapply(ans, paste, collapse = "\r")
        }

        x <- lapply(x, foo, ntip = ntip)

        k <- 0L
        for (i in 1:(n - 1)) {
            y <- x[[i]]
            m1 <- nnode[i]
            for (j in (i + 1):n) {
                z <- x[[j]]
                k <- k + 1L
                res[k] <- m1 + nnode[j] - 2 * sum(z %in% y)
            }
        }
    } else { # method == "score"
        k <- 0L
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                k <- k + 1L
                ## still very slow......
                res[k] <- .dist.topo.score(x[[i]], x[[j]])
            }
        }
    }

    attr(res, "Size") <- n
    attr(res, "Labels") <- nms
    attr(res, "Diag") <- attr(res, "Upper") <- FALSE
    attr(res, "method") <- method
    class(res) <- "dist"
    res
}

.dist.topo.score <- function(x, y)
{
    if (is.null(x$edge.length) || is.null(y$edge.length))
        stop("trees must have branch lengths for branch score distance.")
    nx <- length(x$tip.label)
    x <- reorder.phylo(unroot(x), "postorder")
    y <- reorder.phylo(unroot(y), "postorder")
    ##bp1 <- .Call(bipartition, x$edge, nx, x$Nnode)
    bp1 <- bipartition2(x$edge, nx)
    bp1 <- lapply(bp1, function(xx) sort(x$tip.label[xx]))
    ny <- length(y$tip.label) # fix by Otto Cordero
    ## fix by Tim Wallstrom:
    bp2.tmp <- bipartition2(y$edge, ny)
    ##bp2.tmp <- .Call(bipartition, y$edge, ny, y$Nnode)
    bp2 <- lapply(bp2.tmp, function(xx) sort(y$tip.label[xx]))
    bp2.comp <- lapply(bp2.tmp, function(xx) setdiff(1:ny, xx))
    bp2.comp <- lapply(bp2.comp, function(xx) sort(y$tip.label[xx]))
    ## End
    q1 <- length(bp1)
    q2 <- length(bp2)

    dT <- 0
    found1 <- FALSE
    found2 <- logical(q2)
    found2[1] <- TRUE
    for (i in 2:q1) {
        for (j in 2:q2) {
            if (identical(bp1[[i]], bp2[[j]]) | identical(bp1[[i]], bp2.comp[[j]])) {
                dT <- dT + (x$edge.length[which(x$edge[, 2] == nx + i)] -
                            y$edge.length[which(y$edge[, 2] == ny + j)])^2
                found1 <- found2[j] <- TRUE
                break
            }
        }
        if (found1) found1 <- FALSE
        else dT <- dT + (x$edge.length[which(x$edge[, 2] == nx + i)])^2
    }
    if (!all(found2))
        dT <- dT + sum((y$edge.length[y$edge[, 2] %in% (ny + which(!found2))])^2)
    sqrt(dT)
}

.compressTipLabel <- function(x, ref = NULL)
{
    ## 'x' is a list of objects of class "phylo" possibly with no class
    if (!is.null(attr(x, "TipLabel"))) return(x)
    if (is.null(ref)) ref <- x[[1]]$tip.label
    n <- length(ref)
    if (length(unique(ref)) != n)
        stop("some tip labels are duplicated in tree no. 1")

    ## serious improvement by Joseph W. Brown!
    relabel <- function (y) {
        label <- y$tip.label
        if (!identical(label, ref)) {
            if (length(label) != length(ref))
                stop("one tree has a different number of tips")
            ilab <- match(label, ref)
            if (any(is.na(ilab)))
                stop("one tree has different tip labels")
            ie <- match(1:n, y$edge[, 2])
            y$edge[ie, 2] <- ilab
        }
        y$tip.label <- NULL
        y
    }
    x <- unclass(x) # another killer improvement by Tucson's hackathon (1/2/2013)
    x <- lapply(x, relabel)
    attr(x, "TipLabel") <- ref
    class(x) <- "multiPhylo"
    x
}

prop.part <- function(..., check.labels = TRUE)
{
    obj <- .getTreesFromDotdotdot(...)
    ntree <- length(obj)
    if (ntree == 1) check.labels <- FALSE
    if (check.labels) obj <- .compressTipLabel(obj) # fix by Klaus Schliep (2011-02-21)
    class(obj) <- NULL # fix by Klaus Schliep (2014-03-06)
    for (i in 1:ntree) storage.mode(obj[[i]]$Nnode) <- "integer"
    class(obj) <- "multiPhylo"
    obj <- reorder(obj, "postorder")
# the following line should not be necessary any more
#    obj <- .uncompressTipLabel(obj) # fix a bug (2010-11-18)
    nTips <- length(obj[[1]]$tip.label)
    clades <- prop_part2(obj, nTips)
    attr(clades, "labels") <- obj[[1]]$tip.label
    clades
}

print.prop.part <- function(x, ...)
{
    if (is.null(attr(x, "labels"))) {
        for (i in 1:length(x)) {
            cat("==>", attr(x, "number")[i], "time(s):")
            print(x[[i]], quote = FALSE)
        }
    } else {
        for (i in 1:length(attr(x, "labels")))
          cat(i, ": ", attr(x, "labels")[i], "\n", sep = "")
        cat("\n")
        for (i in 1:length(x)) {
            cat("==>", attr(x, "number")[i], "time(s):")
            print(x[[i]], quote = FALSE)
        }
    }
}

summary.prop.part <- function(object, ...) attr(object, "number")

plot.prop.part <- function(x, barcol = "blue", leftmar = 4, col = "red", ...)
{
    if (is.null(attr(x, "labels")))
      stop("cannot plot this partition object; see ?prop.part for details.")
    L <- length(x)
    n <- length(attr(x, "labels"))
    layout(matrix(1:2, 2, 1), heights = c(1, 3))
    par(mar = c(0.1, leftmar, 0.1, 0.1))
    one2L <- seq_len(L)
    plot(one2L - 0.5, attr(x, "number"), type = "h", col = barcol, xlim = c(0, L),
         xaxs = "i", xlab = "", ylab = "Frequency", xaxt = "n", bty = "n", ...)
    M <- matrix(0L, L, n)
    for (i in one2L) M[i, x[[i]]] <- 1L
    image.default(one2L, 1:n, M, col = c("white", col), xlab = "", ylab = "", yaxt = "n")
    mtext(attr(x, "labels"), side = 2, at = 1:n, las = 1)
}

### by Klaus (2016-03-23):
prop.clades <- function(phy, ..., part = NULL, rooted = FALSE)
{
    if (is.null(part)) {
        obj <- .getTreesFromDotdotdot(...)
        ## avoid double counting of edges if trees are rooted
        if (!rooted) obj <- lapply(obj, unroot)
        part <- prop.part(obj, check.labels = TRUE)
    }
    LABS <- attr(part, "labels")
    if (!identical(phy$tip.label, LABS)) {
        i <- match(phy$tip.label, LABS)
        j <- match(seq_len(Ntip(phy)), phy$edge[, 2])
        phy$edge[j, 2] <- i
        phy$tip.label <- LABS
    }
    bp <- prop.part(phy)
    if (!rooted) {
        ## avoid messing up the order and length if phy is rooted in some cases
        bp <- ONEwise(bp)
        part <- postprocess.prop.part(part)
    }
    pos <- match(bp, part)
    tmp <- which(!is.na(pos))
    n <- rep(NA_real_, phy$Nnode)
    n[tmp] <- attr(part, "number")[pos[tmp]]
    n
}

boot.phylo <-
    function(phy, x, FUN, B = 100, block = 1,
             trees = FALSE, quiet = FALSE,
             rooted = is.rooted(phy), jumble = TRUE,
             mc.cores = 1)
{
    if (is.null(dim(x)) || length(dim(x)) != 2)
        stop("the data 'x' must have two dimensions (e.g., a matrix or a data frame)")

    if (anyDuplicated(rownames(x)))
        stop("some labels are duplicated in the data: you won't be able to analyse tree bipartitions")

    boot.tree <- vector("list", B)
    y <- nc <- ncol(x)
    nr <- nrow(x)

    if (block > 1) {
        a <- seq(1, nc - 1, block)
        b <- seq(block, nc, block)
        y <- mapply(":", a, b, SIMPLIFY = FALSE)
        getBootstrapIndices <- function() unlist(sample(y, replace = TRUE))
    } else getBootstrapIndices <- function() sample.int(y, replace = TRUE)

    if (!quiet) {
        prefix <- "\rRunning bootstraps:      "
        suffix <- paste("/", B)
        updateProgress <- function(i) cat(prefix, i, suffix)
    }

    if (mc.cores == 1) {
        for (i in 1:B) {
            boot.samp <- x[, getBootstrapIndices()]
            if (jumble) boot.samp <- boot.samp[sample.int(nr), ]
            boot.tree[[i]] <- FUN(boot.samp)
            if (!quiet && !(i %% 100)) updateProgress(i)
        }
    } else {
        if (!quiet) cat("Running parallel bootstraps...")
        foo <- function(i) {
            boot.samp <- x[, getBootstrapIndices()]
            if (jumble) boot.samp <- boot.samp[sample.int(nr), ]
            FUN(boot.samp)
        }
        boot.tree <- mclapply(1:B, foo, mc.cores = mc.cores)
        if (!quiet) cat(" done.")
    }

    if (!quiet) cat("\nCalculating bootstrap values...")

    ## sort labels after mixed them up
    if (jumble) {
        boot.tree <- .compressTipLabel(boot.tree, ref = phy$tip.label)
        boot.tree <- .uncompressTipLabel(boot.tree)
        boot.tree <- unclass(boot.tree) # otherwise countBipartitions crashes
    }
    class(boot.tree) <- "multiPhylo"
    if (rooted) {
        pp <- prop.part(boot.tree)
        ans <- prop.clades(phy, part = pp, rooted = rooted)
    } else {
        phy <- reorder(phy, "postorder")
        ints <- phy$edge[, 2] > Ntip(phy)
        ans <- countBipartitions(phy, boot.tree)
        ans <- c(B, ans[order(phy$edge[ints, 2])])
    }

    if (!quiet) cat(" done.\n")

    if (trees) ans <- list(BP = ans, trees = boot.tree)
    ans
}

### The next function transforms an object of class "prop.part" so
### that the vectors which are identical in terms of split are aggregated.
### For instance if n = 5 tips, 1:2 and 3:5 actually represent the same
### split though they are different clades. The aggregation is done
### arbitrarily. The call to ONEwise() insures that all splits include
### the first tip.
### (rewritten by Klaus)
postprocess.prop.part <- function(x)
{
    w <- attr(x, "number")
    labels <- attr(x, "labels")

    x <- ONEwise(x)
    drop <- duplicated(x)

    if (any(drop)) {
        ind1 <- match(x[drop], x)
        ind2 <- which(drop)
        for (i in seq_along(ind2))
            w[ind1[i]] <- w[ind1[i]] + w[ind2[i]]
        x <- x[!drop]
        w <- w[!drop]
    }
    attr(x, "number") <- w
    attr(x, "labels") <- labels
    class(x) <- "prop.part"
    x
}

### This function changes an object of class "prop.part" so that they
### all include the first tip. For instance if n = 5 tips, 3:5 is
### changed to 1:2.
ONEwise <- function(x)
{
    v <- seq_along(x[[1L]])
    for (i in 2:length(x)) {
        y <- x[[i]]
        if (y[1] != 1) x[[i]] <- v[-y]
    }
    x
}

consensus <- function(..., p = 1, check.labels = TRUE)
{
    foo <- function(ic, node) {
        ## ic: index of 'pp'
        ## node: node number in the final tree
        pool <- pp[[ic]]
        if (ic < m) {
            for (j in (ic + 1):m) {
                wh <- match(pp[[j]], pool)
                if (!any(is.na(wh))) {
                    edge[pos, 1] <<- node
                    pool <- pool[-wh]
                    edge[pos, 2] <<- nextnode <<- nextnode + 1L
                    pos <<- pos + 1L
                    foo(j, nextnode)
                }
            }
        }
        size <- length(pool)
        if (size) {
            ind <- pos:(pos + size - 1)
            edge[ind, 1] <<- node
            edge[ind, 2] <<- pool
            pos <<- pos + size
        }
    }

    obj <- .getTreesFromDotdotdot(...)

    if (!is.null(attr(obj, "TipLabel")))
        labels <- attr(obj, "TipLabel")
    else {
        labels <- obj[[1]]$tip.label
        if (check.labels) obj <- .compressTipLabel(obj)
    }
    ntree <- length(obj)
    ## Get all observed partitions and their frequencies:
    pp <- prop.part(obj, check.labels = FALSE)
    ## Drop the partitions whose frequency is less than 'p':
    if (p == 0.5) p <- 0.5000001 # avoid incompatible splits
    pp <- pp[attr(pp, "number") >= p * ntree]
    ## Get the order of the remaining partitions by decreasing size:
    ind <- order(lengths(pp), decreasing = TRUE)
    pp <- pp[ind]
    n <- length(labels)
    m <- length(pp)
    edge <- matrix(0L, n + m - 1, 2)
    if (m == 1) {
        edge[, 1] <- n + 1L
        edge[, 2] <- 1:n
    } else {
        nextnode <- n + 1L
        pos <- 1L
        foo(1, nextnode)
    }
    structure(list(edge = edge, tip.label = labels,
              Nnode = m), class = "phylo")
}
