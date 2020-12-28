## rtt.R (2015-07-16)

##   Root a tree by root-to-tip regression

## Copyright (c) 2014-2015, Rosemary McCloskey, BC Centre for Excellence in HIV/AIDS

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

rtt <- function (t, tip.dates, ncpu = 1, objective = "correlation",
    opt.tol = .Machine$double.eps^0.25)  {
    if (objective == "correlation")
        objective <- function(x, y) cor.test(y, x)$estimate
    else if (objective == "rsquared")
        objective <- function(x, y) summary(lm(y ~ x))$r.squared
    else if (objective == "rms")
        objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
    else stop("objective must be one of \"correlation\", \"rsquared\", or \"rms\"")

    ut <- unroot(t)
    dist <- dist.nodes(ut)[, 1:(ut$Nnode + 2)]

    f <- function (x, parent, child) {
        edge.dist <- x * dist[parent, ] + (1 - x) * dist[child,]
        objective(tip.dates, edge.dist)
    }

    obj.edge <- if (ncpu > 1)
        unlist(mclapply(1:nrow(ut$edge), function (e) {
            opt.fun <- function (x) f(x, ut$edge[e,1], ut$edge[e,2])
            optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
        }, mc.cores=ncpu))
    else apply(ut$edge, 1, function (e) {
        opt.fun <- function (x) f(x, e[1], e[2])
        optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
    })

    best.edge <- which.max(obj.edge)

    best.edge.parent <- ut$edge[best.edge, 1]
    best.edge.child <- ut$edge[best.edge, 2]
    best.edge.length <- ut$edge.length[best.edge]

    opt.fun <- function (x) f(x, best.edge.parent, best.edge.child)
    best.pos <- optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum

    new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root",
        edge.length = 1, Nnode = 1L, root.edge = 1)
    class(new.root) <- "phylo"
    ut <- bind.tree(ut, new.root, where = best.edge.child, position = best.pos *
        best.edge.length)
    ut <- collapse.singles(ut)
    ut <- root(ut, "new.root")
    drop.tip(ut, "new.root")
}
