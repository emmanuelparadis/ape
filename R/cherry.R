## cherry.R (2009-05-10)

##     Number of Cherries and Null Models of Trees

## Copyright 2002-2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

cherry <- function(phy)
{
    if (!inherits(phy, "phylo")) stop('object "phy" is not of class "phylo"')
    n <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != n - 1) stop('"phy" is not fully dichotomous')
    if (n < 4) stop("not enough tips in your phylogeny for this analysis")
    cherry <- sum(tabulate(phy$edge[, 1][phy$edge[, 2] <= n]) == 2)
    small.n <- n < 20
    if (small.n) {
        P.yule <- f.cherry.yule(n, cherry)
        P.uniform <- f.cherry.uniform(n, cherry)
    }
    else {
        P.yule <- 2*(1 - pnorm(abs(cherry - n/3)/sqrt(2*n/45)))
        mu.unif <- n*(n - 1)/(2*(2*n - 5))
        sigma2.unif <- n*(n - 1)*(n - 4)*(n - 5)/(2*(2*n - 5)^2 * (2*n -7))
        P.uniform <- 2*(1 - pnorm(abs(cherry - mu.unif)/sqrt(sigma2.unif)))
    }
    cat("\nAnalysis of the Number of Cherries in a Tree\n\n")
    cat("Phylogenetic tree:", deparse(substitute(phy)), "\n")
    cat("Number of tips:", n, "\n")
    cat("Number of cherries:", cherry, "\n\n")
    cat("Null hypothesis: Yule model\n")
    cat("    P-value =", round(P.yule, 4), "\n\n")
    cat("Null hypothesis: uniform model\n")
    cat("    P-value =", round(P.uniform, 4), "\n\n")
    if (!small.n) cat("(P-values were computed using normal approximations)\n")
}

f.cherry.yule <- function(n, k)
{
    if (k == 0 || k > floor(n/2)) 0 else if (n == 4) if (k == 1) 2/3 else if (k == 2) 1/3 else 0
    else (1 - 2*(k - 1)/(n - 1))*f.cherry.yule(n - 1, k - 1) +
        2*k/(n - 1)*f.cherry.yule(n - 1, k)
}

f.cherry.uniform <- function(n, k)
{
    if (k == 0 || k > floor(n/2)) 0 else if (n == 4) if (k == 1) 4/5 else if (k == 2) 1/5 else 0
    else if (k == 1) 0 else (gamma(n + 1)*gamma(n - 2 + 1)*gamma(n - 4 + 1) * 2^(n-2*k)) /
        (gamma(n - 2*k + 1)*gamma(2*n - 4 + 1)*gamma(k + 1)*gamma(k - 2 + 1))
}
