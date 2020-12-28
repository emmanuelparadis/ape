## plot.popsize.R (2004-07-4) modified by EP (2019-01-29)

##   Plot population size in dependence of time

## Copyright 2004 Rainer Opgen-Rhein and Korbinian Strimmer

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

plot.popsize <-
    function(x, show.median = TRUE, show.years = FALSE,
             subst.rate, present.year, xlab = NULL,
             ylab = "Effective population size",
             log = "y", ...)
{
    ylim <- range(x[, 2:5], na.rm = TRUE)
    x1 <- x[, 1]
    if (show.years) {
        x1 <- -x1/subst.rate + present.year
        if (is.null(xlab)) xlab <- "Time (years)"
    } else {
        if (is.null(xlab))
            xlab <- "Time (past to present in units of substitutions)"
    }
    xlim <- range(x1, na.rm = TRUE)

    j <- if (show.median) 3 else 2

    plot(x1, x[, j], type = "s", xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, log = log, lwd = 2.5, ...)

    lines(x1, x[, 4], ...)
    lines(x1, x[, 5], ...)
}

lines.popsize <- function(x, show.median = TRUE, show.years = FALSE,
                          subst.rate, present.year, ...)
{
    x1 <- x[, 1]
    if (show.years) x1 <- -x1/subst.rate + present.year
    j <- if (show.median) 3 else 2
    lines(x1, x[, j], lwd = 2.5, ...)
    lines(x1, x[, 4], ...)
    lines(x1, x[, 5], ...)
}
