## compar.ou.R (2010-11-04)

##   Ornstein--Uhlenbeck Model for Continuous Characters

## Copyright 2005-2010 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

compar.ou <- function(x, phy, node = NULL, alpha = NULL)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo".')
    if (!is.numeric(x)) stop("'x' must be numeric.")
    if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label)) x <- x[phy$tip.label]
        else warning('the names of argument "x" and the tip labels of the tree did not match: the former were ignored in the analysis.')
    }
    n <- length(phy$tip.label)
    root <- n + 1L
    if (is.null(node)) node <- numeric(0)
    if (is.character(node)) {
        if (is.null(phy$node.label))
            stop("argument 'node' is character but 'phy' has no node labels")
        node <- match(node, phy$node.label) + n
        phy$node.label <- NULL
    }
    if (root %in% node) node <- node[-1]
    bt <- branching.times(phy)
    Tmax <- bt[1]
    Wend <- matrix(0, n, length(node) + 1)
    colnames(Wend) <- c(names(sort(bt[node - n])), as.character(root))
    Wstart <- Wend
    Wstart[, ncol(Wstart)] <- Tmax
    root2tip <- .Call(seq_root2tip, phy$edge, n, phy$Nnode)
    for (i in 1:n) {
        last.change <- names(Tmax)
        for (j in root2tip[[i]]) {
            if (j %in% node) {
                jb <- as.character(j)
                Wend[i, last.change] <- Wstart[i, jb] <- bt[jb]
                last.change <- jb
            }
        }
    }
    W <- cophenetic.phylo(phy)
    dev <- function(p) {
        alpha <- p[1]
        sigma2 <- p[2]
        if (sigma2 <= 0) return(1e100)
        theta <- p[-(1:2)]
        ## fixed a bug below: must be '%*% theta' instead of '* theta' (2010-03-15)
        M <- rowSums((exp(-alpha * Wend) - exp(-alpha * Wstart)) %*% theta)
        V <- exp(-alpha * W) * (1 - exp(-2 * alpha * (Tmax - W/2)))
        R <- chol(V) # correction by Cecile Ane (2010-11-04)
        n * log(2 * pi * sigma2) + 2 * sum(log(diag(R))) +
            (t(x - M) %*% chol2inv(R) %*% (x - M)) / sigma2
    }

    out <- if (is.null(alpha))
        nlm(function(p) dev(p),
            p = c(0.1, 1, rep(mean(x), ncol(Wstart))), hessian = TRUE)
    else nlm(function(p) dev(c(alpha, p)),
             p = c(1, rep(mean(x), ncol(Wstart))), hessian = TRUE)

    ## if alpha is estimated it may be that the Hessian matrix has the
    ## corresponding column and row filled with 0, making solve() fail
    se <- if (is.null(alpha) && all(out$hessian[1, ] == 0))
        c(NA, sqrt(diag(solve(out$hessian[-1, -1]))))
    else sqrt(diag(solve(out$hessian)))

    para <- cbind(out$estimate, se)
    nms <- c("sigma2", paste("theta", 1:ncol(Wstart), sep = ""))
    if (is.null(alpha)) nms <- c("alpha", nms)
    dimnames(para) <- list(nms, c("estimate", "stderr"))
    obj <- list(deviance = out$minimum, para = para, call = match.call())
    class(obj) <- "compar.ou"
    obj
}
