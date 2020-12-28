## alex.R (2017-04-18)

##   Alignment Explorer With Multiple Devices

## Copyright 2012-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

alex <- function(x, ...)
{
    n <- nrow(x)
    s <- ncol(x)
    devmain <- dev.cur()
    on.exit(dev.set(devmain))
    NEW <- TRUE
    cat("Click on two opposite corners of the zone you want to zoom-in.
Right-click to exit.\n")
    repeat {
        xy <- locator(2)
        if (is.null(xy)) break
        xy$y <- n - xy$y + 1
        xy <- lapply(xy, function(x) sort(round(x)))
        i1 <- xy$y[1L]; i2 <- xy$y[2L]
        j1 <- xy$x[1L]; j2 <- xy$x[2L]
        if (i1 > n || j1 > s) cat("Try again!\n") else {
            if (i1 <= 0) i1 <- 1L
            if (j1 <= 0) j1 <- 1L
            if (i2 > n) i2 <- n
            if (j2 > s) j2 <- s
            if (NEW) {
                dev.new()
                devsub <- dev.cur()
                NEW <- FALSE
            } else dev.set(devsub)
            image(x[i1:i2, j1:j2], xaxt = "n", ...)
            atx <- axTicks(1)
            axis(1, atx, labels = (j1:j2)[atx])
            title(sub = paste("From", sQuote(deparse(substitute(x)))))
            dev.set(devmain)
        }
    }
}
