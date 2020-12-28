## scales.R (2014-08-21)

##   Add a Scale Bar or Axis to a Phylogeny Plot

## add.scale.bar: add a scale bar to a phylogeny plot
## axisPhylo: add a scale axis on the side of a phylogeny plot

## Copyright 2002-2014 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

add.scale.bar <- function(x, y, length = NULL, ask = FALSE,
                          lwd = 1, lcol = "black", ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    direc <- lastPP$direction
    if (is.null(length)) {
        nb.digit <-
          if (direc %in% c("rightwards", "leftwards")) diff(range(lastPP$xx))
          else diff(range(lastPP$yy))
        length <- pretty(c(0, nb.digit) / 6, 1)[2] # by Klaus
    }

    if (ask) {
        cat("\nClick where you want to draw the bar\n")
        x <- unlist(locator(1))
        y <- x[2]
        x <- x[1]
    } else if (missing(x) || missing(y)) {
        if (lastPP$type %in% c("phylogram", "cladogram")) {
            switch(direc,
                   "rightwards" = {
                       x <- 0
                       y <- 1
                   },
                   "leftwards" = {
                       x <- max(lastPP$xx)
                       y <- 1
                   },
                   "upwards" = {
                       x <- max(lastPP$xx)
                       y <- 0
                   },
                   "downwards" = {
                       x <- 1
                       y <- max(lastPP$yy)
                   })
        } else {
            direc <- "rightwards" # just to be sure for below
            x <- lastPP$x.lim[1]
            y <- lastPP$y.lim[1]
        }
    }

    switch(direc,
           "rightwards" = {
               segments(x, y, x + length, y, col = lcol, lwd = lwd)
               text(x + length * 1.1, y, as.character(length), adj = c(0, 0.5), ...)
           },
           "leftwards" = {
               segments(x - length, y, x, y, col = lcol, lwd = lwd)
               text(x - length * 1.1, y, as.character(length), adj = c(1, 0.5), ...)
           },
           "upwards" = {
               segments(x, y, x, y + length, col = lcol, lwd = lwd)
               text(x, y + length * 1.1, as.character(length), adj = c(0, 0.5), srt = 90, ...)
           },
           "downwards" = {
               segments(x, y - length, x, y, col = lcol, lwd = lwd)
               text(x, y - length * 1.1, as.character(length), adj = c(0, 0.5), srt = 270, ...)
           })
}

axisPhylo <- function(side = 1, root.time = NULL, backward = TRUE, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    type <- lastPP$type

    if (type == "unrooted")
        stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
    if (type == "radial")
        stop("axisPhylo() not meaningful for this type of plot")

    if (is.null(root.time)) root.time <- lastPP$root.time

    if (type %in% c("phylogram", "cladogram")) {
        xscale <-
            if (lastPP$direction %in% c("rightwards", "leftwards")) range(lastPP$xx)
            else range(lastPP$yy)

        tmp <- lastPP$direction %in% c("leftwards", "downwards")

        tscale <- c(0, xscale[2] - xscale[1])
        if (xor(backward, tmp)) tscale <- tscale[2:1]
        if (!is.null(root.time)) {
            tscale <- tscale + root.time
            if (backward) tscale <- tscale - xscale[2]
        }

        ## the linear transformation between the x-scale and the time-scale:
        beta <- diff(xscale) / diff(tscale)
        alpha <- xscale[1] - beta * tscale[1]

        lab <- pretty(tscale)
        x <- beta * lab + alpha
        axis(side = side, at = x, labels = lab, ...)
    } else { # type == "fan"
        n <- lastPP$Ntip
        xx <- lastPP$xx[1:n]; yy <- lastPP$yy[1:n]
        r0 <- max(sqrt(xx^2 + yy^2))
        firstandlast <- c(1, n)
        theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
        x0 <- r0 * cos(theta0); y0 <- r0 * sin(theta0)
        inc <- diff(pretty(c(0, r0))[1:2])
        srt <- 360*theta0/(2*pi)
        coef <- -1
        if (abs(srt) > 90) {
            srt <- srt + 180
            coef <- 1
        }
        len <- 0.025 * r0 # the length of tick marks
        r <- r0
        while (r > 1e-8) {
            x <- r * cos(theta0); y <- r * sin(theta0)
            if (len/r < 1) {
                ra <- sqrt(len^2 + r^2); thetaa <- theta0 + coef * asin(len/r)
                xa <- ra * cos(thetaa); ya <- ra * sin(thetaa)
                segments(xa, ya, x, y)
                text(xa, ya, r0 - r, srt = srt, adj = c(0.5, 1.1), ...)
            }
            r <- r - inc
        }
        segments(x, y, x0, y0)
    }
}
