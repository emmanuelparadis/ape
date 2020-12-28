## is.compatible.R (2017-06-03)

##   Check Compatibility of Splits

## Copyright 2011 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

is.compatible <- function(obj) UseMethod("is.compatible")

is.compatible.bitsplits <- function(obj)
{
    m <- obj$matsplit
    n <- ncol(m)
    ntaxa <- length(obj$labels)
    for (i in 1:(n - 1))
        for (j in (i + 1):n)
            if (!arecompatible(m[, i], m[, j], ntaxa))
                return(FALSE)
    TRUE
}

arecompatible <-function(x, y, n)
{
    msk <- !as.raw(2^(8 - (n %% 8)) - 1)

    foo <- function(v) {
        lv <- length(v)
        v[lv] <- v[lv] & msk
        as.integer(all(v == as.raw(0)))
    }

    nE <- foo(x & y) + foo(x & !y) + foo(!x & y) + foo(!x & !y)
    if (nE >= 1) TRUE else FALSE
}
