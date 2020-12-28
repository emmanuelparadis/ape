## birthdeath.R (2012-04-20)

##   Estimation of Speciation and Extinction Rates
##             with Birth-Death Models

## birthdeath: standard model
## bd.ext: extended version

## Copyright 2002-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

birthdeath <- function(phy)
{
    if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo"')
    N <- length(phy$tip.label)
    x <- c(NA, branching.times(phy))
    dev <- function(a, r) {
        if (r < 0 || a > 1) return(1e100)
        -2 * (lfactorial(N - 1)
              + (N - 2) * log(r)
              + r * sum(x[3:N])
              + N * log(1 - a)
              - 2 * sum(log(exp(r * x[2:N]) - a)))
    }
    out <- nlm(function(p) dev(p[1], p[2]), c(0.1, 0.2), hessian = TRUE)
    if (out$estimate[1] < 0) {
        out <- nlm(function(p) dev(0, p), 0.2, hessian = TRUE)
        para <- c(0, out$estimate)
        inv.hessian <- try(solve(out$hessian))
        se <-
            if (class(inv.hessian) == "try-error") NA
            else sqrt(diag(inv.hessian))
        se <- c(0, se)
    }
    else {
        para <- out$estimate
        inv.hessian <- try(solve(out$hessian))
        se <-
            if (class(inv.hessian) == "try-error") c(NA, NA)
            else sqrt(diag(inv.hessian))
    }
    Dev <- out$minimum

    ## 95% profile likelihood CIs

    ## which: index of the parameter (1 or 2)
    ## s: sign of the increment (-1 or +1)
    foo <- function(which, s) {
        i <- 0.1

        if (which == 1) {
            p <- para[1] + s * i
            bar <- function() dev(p, para[2])
        } else { # which == 2
            p <- para[2] + s * i
            bar <- function() dev(para[1], p)
        }

        while (i > 1e-9) {
            while (bar() < Dev + 3.84) p <- p + s * i
            p <- p - s * i
            i <- i / 10
        }
        p
    }

    CI <- mapply(foo, c(1, 2, 1, 2), c(-1, -1, 1, 1))
    dim(CI) <- c(2, 2)

    names(para) <- names(se) <- rownames(CI) <- c("d/b", "b-d")
    colnames(CI) <- c("lo", "up")
    obj <- list(tree = deparse(substitute(phy)), N = N,
                dev = Dev, para = para, se = se, CI = CI)
    class(obj) <- "birthdeath"
    obj
}

print.birthdeath <- function(x, ...)
{
    cat("\nEstimation of Speciation and Extinction Rates\n")
    cat("            with Birth-Death Models\n\n")
    cat("     Phylogenetic tree:", x$tree, "\n")
    cat("        Number of tips:", x$N, "\n")
    cat("              Deviance:", x$dev, "\n")
    cat("        Log-likelihood:", -(x$dev)/2, "\n")
    cat("   Parameter estimates:\n")
    cat("      d / b =", x$para[1], "  StdErr =", x$se[1], "\n")
    cat("      b - d =", x$para[2], "  StdErr =", x$se[2], "\n")
    cat("   (b: speciation rate, d: extinction rate)\n")
    cat("   Profile likelihood 95% confidence intervals:\n")
    cat("      d / b: [", x$CI[1, 1], ", ", x$CI[1, 2], "]", "\n", sep = "")
    cat("      b - d: [", x$CI[2, 1], ", ", x$CI[2, 2], "]", "\n\n", sep = "")
}

bd.ext <- function(phy, S, conditional = TRUE)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    if (!is.null(names(S))) {
        if (all(names(S) %in% phy$tip.label)) S <- S[phy$tip.label]
        else warning('the names of argument "S" and the tip labels
did not match: the former were ignored.')
    }
    N <- length(S)
    x <- branching.times(phy)
    x <- c(x[1], x)
    trm.br <- phy$edge.length[phy$edge[, 2] <= N]
    if (conditional) {
        dev <- function(a, r) {
            if (a >= 1 || a < 0 || r <= 0) return(1e50)
            ert <- exp(r * trm.br)
            zeta <- (ert - 1)/(ert - a)
            -2 * (lfactorial(N - 1)
                  + (N - 2) * log(r)
                  + N * log(1 - a)
                  + 2 * r * sum(x[2:N])
                  - 2 * sum(log(exp(r * x[2:N]) - a))
                  + sum(log(1 - zeta) + (S - 1)*log(zeta)))
        }
    } else {
        dev <- function(a, r) {
            if (a >= 1 || a < 0 || r <= 0) return(1e50)
            -2 * (lfactorial(N - 1)
                  + (N - 2) * log(r)
                  + (3 * N) * log(1 - a)
                  + 2 * r * sum(x[2:N])
                  - 2 * sum(log(exp(r * x[2:N]) - a))
                  + r * sum(trm.br)
                  + sum((S - 1) * log(exp(r * trm.br) - 1))
                  - sum((S + 1) * log(exp(r * trm.br) - a)))
        }
    }
    out <- nlm(function(p) dev(p[1], p[2]), c(0.1, 0.2), hessian = TRUE)
    para <- out$estimate
    se <- sqrt(diag(solve(out$hessian)))
    Dev <- out$minimum
    cat("\nExtended Version of the Birth-Death Models to\n")
    cat("    Estimate Speciation and Extinction Rates\n\n")
    cat("    Data: phylogenetic:", deparse(substitute(phy)), "\n")
    cat("             taxonomic:", deparse(substitute(S)), "\n")
    cat("        Number of tips:", N, "\n")
    cat("              Deviance:", Dev, "\n")
    cat("        Log-likelihood:", -Dev/2, "\n")
    cat("   Parameter estimates:\n")
    cat("      d / b =", para[1], "  StdErr =", se[1], "\n")
    cat("      b - d =", para[2], "  StdErr =", se[2], "\n")
    cat("   (b: speciation rate, d: extinction rate)\n")
}
