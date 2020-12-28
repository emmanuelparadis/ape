## rTrait.R (2014-03-06)

##   Trait Evolution

## Copyright 2010-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

rTraitDisc <-
    function(phy, model = "ER", k = if (is.matrix(model)) ncol(model) else 2,
             rate = 0.1, states = LETTERS[1:k], freq = rep(1/k, k),
             ancestor = FALSE, root.value = 1, ...)
{
    if (is.null(phy$edge.length))
        stop("tree has no branch length")
    if (any(phy$edge.length < 0))
        stop("at least one branch length negative")

    if (is.character(model)) {
        switch(toupper(model), "ER" = {
                   if (length(rate) != 1)
                       stop("`rate' must have one element")
                   Q <- matrix(rate, k, k)
               }, "ARD" = {
                   if (length(rate) != k*(k - 1))
                       stop("`rate' must have k(k - 1) elements")
                   Q <- matrix(0, k, k)
                   Q[col(Q) != row(Q)] <- rate
               }, "SYM" = {
                   if (length(rate) != k*(k - 1)/2)
                       stop("`rate' must have k(k - 1)/2 elements")
                   Q <- matrix(0, k, k)
                   sel <- col(Q) < row(Q)
                   Q[sel] <- rate
                   Q <- t(Q)
                   Q[sel] <- rate
               })
    }
    if (is.matrix(model)) {
        Q <- model
        if (ncol(Q) != nrow(Q))
            stop("the matrix given as `model' must be square")
    }

    phy <- reorder(phy, "postorder")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    x <- integer(n + phy$Nnode)
    x[ROOT] <- as.integer(root.value)

    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    el <- phy$edge.length

    if (is.function(model)) {
        environment(model) <- environment() # to find 'k'
        for (i in N:1) x[des[i]] <- model(x[anc[i]], el[i], ...)
    } else {
        freq <- rep(freq, each = k)
        Q <- Q * freq
        diag(Q) <- 0
        diag(Q) <- -rowSums(Q)
        for (i in N:1) {
            p <- matexpo(Q * el[i])[x[anc[i]], ]
            x[des[i]] <- sample.int(k, size = 1, FALSE, prob = p)
        }
    }

    if (ancestor) {
        if (is.null(phy$node.label)) phy <- makeNodeLabel(phy)
        names(x) <- c(phy$tip.label, phy$node.label)
    } else {
        x <- x[1:n]
        names(x) <- phy$tip.label
    }
    class(x) <- "factor"
    levels(x) <- states
    x
}

rTraitCont <-
    function(phy, model = "BM", sigma = 0.1, alpha = 1, theta = 0,
             ancestor = FALSE, root.value = 0, ...)
{
    if (is.null(phy$edge.length))
        stop("tree has no branch length")
    if (any(phy$edge.length < 0))
        stop("at least one branch length negative")

    phy <- reorder(phy, "postorder")
    n <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L
    x <- numeric(n + phy$Nnode)
    x[ROOT] <- root.value

    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    el <- phy$edge.length

    if (is.function(model)) {
        environment(model) <- environment()
        for (i in N:1) x[des[i]] <- model(x[anc[i]], el[i], ...)
    } else {
        model <- pmatch(toupper(model), c("BM", "OU"))
        if (length(sigma) == 1) sigma <- rep(sigma, N)
        else if (length(sigma) != N)
            stop("'sigma' must have one or Nedge(phy) elements")
        if (model == 2) { # "OU"
            if (length(alpha) == 1) alpha <- rep(alpha, N)
            else if (length(alpha) != N)
                stop("'alpha' must have one or Nedge(phy) elements")
            if (length(theta) == 1) theta <- rep(theta, N)
            else if (length(theta) != N)
                stop("'theta' must have one or Nedge(phy) elements")
        }
        x <- .C(C_rTraitCont, as.integer(model), as.integer(N),
                as.integer(anc - 1L), as.integer(des - 1L), el,
                as.double(sigma), as.double(alpha), as.double(theta),
                x = x, NAOK = TRUE)$x
    }

    if (ancestor) {
        if (is.null(phy$node.label)) phy <- makeNodeLabel(phy)
        names(x) <- c(phy$tip.label, phy$node.label)
    } else {
        x <- x[1:n]
        names(x) <- phy$tip.label
    }
    x
}

rTraitMult <-
    function(phy, model, p = 1, root.value = rep(0, p), ancestor = FALSE,
             asFactor = NULL, trait.labels = paste("x", 1:p, sep = ""), ...)
{
    phy <- reorder(phy, "postorder")
    n <- length(phy$tip.label)
    m <- phy$Nnode
    N <- dim(phy$edge)[1]
    ROOT <- n + 1L

    x <- matrix(0, n + m, p)
    x[ROOT, ] <- root.value

    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]

    el <- phy$edge.length
    if (is.null(el)) el <- numeric(N)

    environment(model) <- environment() # to find 'p'

    for (i in N:1) x[des[i], ] <- model(x[anc[i], ], el[i], ...)

    if (ancestor) {
        if (is.null(phy$node.label)) phy <- makeNodeLabel(phy)
        rownames(x) <- c(phy$tip.label, phy$node.label)
    } else {
        x <- x[1:n, , drop = FALSE]
        rownames(x) <- phy$tip.label
    }
    x <- as.data.frame(x)
    names(x) <- trait.labels
    if (!is.null(asFactor)) {
        for (i in asFactor) {
            y <- x[, i]
            x[, i] <- factor(y, labels = LETTERS[1:length(unique(y))])
        }
    }
    x
}
