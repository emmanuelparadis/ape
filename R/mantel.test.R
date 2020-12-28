## mantel.test.R (2019-02-25)

##   Mantel Test for Similarity of Two Matrices

## Copyright 2002-2011 Ben Bolker and Julien Claude, 2019 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

perm.rowscols <- function(m1, n)
{
    s <- sample(1:n)
    m1[s, s]
}

## calculate the Mantel z-statistic for two square matrices m1 and m2
## old code:
## mant.zstat <- function(m1, m2) sum(lower.triang(m1 * m2))
## modified by EP following suggestion by Andrzej Galecki (2018-02-07)
mant.zstat <- function(m1, m2) {
    diag(m1) <- diag(m2) <- 0 # in case the diagonals are not 0
    sum(m1 * m2)/2
}

mantel.test <- function(m1, m2, nperm = 999, graph = FALSE,
                        alternative = "two.sided", ...)
{
    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    n <- nrow(m1)
    realz <- mant.zstat(m1, m2)
    nullstats <- replicate(nperm, mant.zstat(m1, perm.rowscols(m2, n)))
    pval <- switch(alternative,
                   "two.sided" = 2 * min(sum(nullstats >= realz), sum(nullstats <= realz)),
                   "less" = sum(nullstats <= realz),
                   "greater" = sum(nullstats >= realz))
    pval <- (pval + 1) / (nperm + 1) # 'realz' is included in 'nullstats'
    if (alternative == "two.sided" && pval > 1) pval <- 1
    if (graph) {
        plot(density(nullstats), type = "l", ...)
        abline(v = realz)
    }
    list(z.stat = realz, p = pval, alternative = alternative)
}
