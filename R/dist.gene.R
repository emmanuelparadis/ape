## dist.gene.R (2012-04-02)

##   Pairwise Distances from Genetic Data

## Copyright 2002-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dist.gene <-
    function(x, method = "pairwise", pairwise.deletion = FALSE,
             variance = FALSE)
{
    if (is.data.frame(x)) x <- as.matrix(x) else { # suggestion by Markus Schlegel
        if (!is.matrix(x))
            stop("'x' should be a matrix or a data.frame")
    }
    method <- match.arg(method, c("pairwise", "percentage"))

    if (!pairwise.deletion) {
        ## delete the columns with at least one NA:
        del <- apply(x, 2, function(xx) any(is.na(xx)))
        x <- x[, !del]
    }
    n <- dim(x)
    L <- n[2]
    n <- n[1]
    D <- double(n * (n - 1)/2)
    if (pairwise.deletion) L <- D
    k <- 1L
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            y <- x[i, ] != x[j, ]
            if (pairwise.deletion) L[k] <- sum(!is.na(y))
            D[k] <-  sum(y, na.rm = TRUE)
            k <- k + 1L
        }
    }
    ## L is either a single integer value if pairwise.deletion = FALSE,
    ## or a vector of integers if pairwise.deletion = TRUE

    if (method == "percentage") D <- D/L

    attr(D, "Size") <- n
    attr(D, "Labels") <-  dimnames(x)[[1]]
    attr(D, "Diag") <- attr(D, "Upper") <- FALSE
    attr(D, "call") <- match.call()
    attr(D, "method") <- method
    class(D) <- "dist"

    if (variance) {
        y <- if (method == "pairwise") L else 1
        attr(D, "variance") <- D * (y - D)/L
    }
    D
}
