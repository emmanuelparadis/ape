## delta.plot.R (2010-01-12)

##   Delta Plots

## Copyright 2010 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

delta.plot <- function(X, k = 20, plot = TRUE, which = 1:2)
{
    if (is.matrix(X)) X <- as.dist(X)
    n <- attr(X, "Size")
    if (n < 4) stop("need at least 4 observations")

    ## add a category for the cases delta = 1
    ans <- .C(delta_plot, as.double(X), as.integer(n),
              as.integer(k), integer(k + 1), double(n),
              NAOK = TRUE)
    counts <- ans[[4]]
    ## add the counts of delta=1 to the last category:
    counts[k] <- counts[k] + counts[k + 1]
    counts <- counts[-(k + 1)]

    delta.bar <- ans[[5]]/choose(n - 1, 3)

    if (plot) {
        if (length(which) == 2) layout(matrix(1:2, 1, 2))
        if (1 %in% which) {
            barplot(counts, space = 0, xlab = expression(delta[q]))
            a <- axTicks(1)
            axis(1, at = a, labels = a/k)
        }
        if (2 %in% which)
            plot(delta.bar, type = "h", ylab = expression(bar(delta)))
    }

    invisible(list(counts = counts, delta.bar = delta.bar))
}
