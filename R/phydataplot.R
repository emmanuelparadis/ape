## phydataplot.R (2017-10-04)

##   Annotate Phylogenies

## Copyright 2014-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

polar2rect <- function(r, angle)
    list(x = r * cos(angle), y = r * sin(angle))

rect2polar <- function(x, y)
    list(r = sqrt(x^2 + y^2), angle = atan2(y, x))

.matchDataPhylo <- function(x, phy)
{
    msg <- "'x' has no (row)names: data are assumed to be in the same order than the tips of the tree"
    labs <- phy$tip.label
    if (is.vector(x)) { # also for lists
        if (is.null(names(x))) warning(msg) else x <- x[labs]
    } else {
        if (is.null(rownames(x))) warning(msg) else x <- x[labs, ]
    }
    x
}

ring <- function(x, phy, style = "ring", offset = 1, ...)
{
    style <- match.arg(style, c("ring", "segments", "arrows"))
    x <- .matchDataPhylo(x, phy)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    n <- lastPP$Ntip
    one2n <- seq_len(n)
    tmp <- rect2polar(lastPP$xx[one2n], lastPP$yy[one2n])
    theta <- tmp$angle
    r0 <- max(tmp$r) + offset
    r1 <- r0 + x
    s0 <- polar2rect(rep.int(r0, 100L), seq(0, 2*pi, length.out = 100L))
    s1 <- polar2rect(r1, theta)

    switch(style, ring = {
        if (length(x) < n) x <- rep_len(x, n)
        dx <- dim(x)
        if (is.null(dx)) dim(x) <- dx <- c(n, 1L)
        nc <- dx[2]
        col <- list(...)$col
        if (is.null(col)) col <- "grey"
        if (nc == 1) {
            col <- rep_len(col, n)
        } else {
            colvar <- col
            col <- rep(col[1], n)
        }
        iangle <- min(diff(sort(theta)))
        iangle2 <- iangle / 2
        for (i in one2n) {
            R <- rep(r0, 100)
            THETA <- seq(theta[i] - iangle2, theta[i] + iangle2, length.out = 100)
            xy1 <- polar2rect(R, THETA)
            xy2 <- polar2rect(R + x[i, 1], THETA)
            polygon(c(xy1$x, rev(xy2$x)), c(xy1$y, rev(xy2$y)), col = col[i], border = NA)
            if (nc > 1) {
                for (j in 2:nc) {
                    xy1 <- xy2
                    xy2 <- polar2rect(R + sum(x[i, 1:j]), THETA)
                    polygon(c(xy1$x, rev(xy2$x)), c(xy1$y, rev(xy2$y)), col = colvar[j], border = NA)
                }
            }
        }
        ##polypath(c(s0$x, NA, s0$x), c(s0$y, NA, s1$y), rule = "evenodd",
        ##         border = 1, col = "transparent")
    }, segments = {
        s0 <- polar2rect(rep.int(r0, n), theta)
        segments(s0$x, s0$y, s1$x, s1$y, ...)
    },  arrows = {
        s0 <- polar2rect(rep.int(r0, n), theta)
        fancyarrows(s0$x, s0$y, s1$x, s1$y, ...)
    })
}

phydataplot <- function(x, phy, style = "bars", offset = 1, scaling = 1,
                        continuous = FALSE, width = NULL, legend = "below",
                        funcol = rainbow, ...)
{
    style <- match.arg(style, c("bars", "segments", "image", "arrows", "boxplot", "dotchart", "mosaic"))
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    circular <- if (lastPP$type %in% c("radial", "fan")) TRUE else FALSE

    n <- length(phy$tip.label)
    one2n <- seq_len(n)
    x <- .matchDataPhylo(x, phy)

    if (scaling != 1)
        x <- if (is.list(x)) lapply(x, "*", scaling) else scaling * x

    if (!circular) {
        if (lastPP$direction != "rightwards")
            stop("for the moment, only rightwards trees are supported")
        x0 <- max(lastPP$xx[one2n]) + offset
        if (style %in% c("bars", "segments", "arrows")) x1 <- x0 + x
        y1 <- lastPP$yy[one2n]
        if (style %in% c("bars", "image", "boxplot", "dotchart", "mosaic")) {
            o <- order(y1)
            x <- if (style == "image") x[o, o] else
            if (is.vector(x)) x[o] else x[o, ]
        }
    } else {
        if (style %in% c("image", "boxplot", "dotchart", "mosaic"))
            stop(paste(dQuote(style), "not implemented with circular trees"))
    }

    switch(style, bars = {
        if (circular)
            stop("style = \"bars\" not implemented with circular trees; see function 'ring'")
        if (!is.null(dim(x))) x <- t(x)
        barplot(x, width = 1, add = TRUE, horiz = TRUE, offset = x0,
                axes = FALSE, axisnames = FALSE, space = c(0.5, rep(0, n - 1)), ...)
        px <- pretty(c(0, x))
        axis(1, px + x0, labels = px / scaling, line = 1)
    }, segments = {
        if (circular) ring(x, phy, style, offset, ...)
        else segments(x0, y1, x1, y1, ...)
    }, image = {
        if (inherits(x, "DNAbin"))
            stop('object of class "DNAbin" not supported: use type="mosaic"')
        x1 <- seq(x0, lastPP$x.lim[2], length.out = n)
        image(x1, y1[o], x, add = TRUE, ...)
        mtext(phy$tip.label[o], 1, 1, at = x1, font = lastPP$font,
              cex = lastPP$cex, col = lastPP$tip.color)
    }, arrows = {
        if (circular) ring(x, phy, style, offset, ...)
        else fancyarrows(rep(x0, length(y1)), y1, x1, y1, ...)
    }, boxplot = {
        if (is.matrix(x)) x <- t(x)
        o <- boxplot(x, plot = FALSE)
        mini <- min(o$stats)
        maxi <- max(o$stats)
        if (length(o$out)) { # in case there is no outlier
            mini <- min(o$out, mini)
            maxi <- max(o$out, maxi)
        }
        px <- pretty(c(mini, maxi))
        x0 <- x0 - mini
        o$stats <- o$stats + x0
        o$out <- o$out + x0
        bxp(o, horizontal = TRUE, add = TRUE, axes = FALSE, ...)
        axis(1, px + x0, labels = px / scaling, line = 1)

    }, dotchart = {
        mini <- min(x)
        maxi <- max(x)
        x0 <- x0 - mini
        segments(mini + x0, one2n, maxi + x0, one2n, lty = 3, col = "gray")
        points(x + x0, 1:n, ...)
        px <- pretty(x)
        axis(1, px + x0, labels = px / scaling, line = 1)
    }, mosaic = {
        p <- ncol(x)
        if (is.null(p)) p <- 1L
        if (is.null(width)) {
            x1 <- lastPP$x.lim[2]
            width <- (x1 - x0)/p
        } else x1 <- x0 + width * p
        xx <- seq(x0, x1, width)
        xl <- rep(xx[-length(xx)], each = n)
        yb <- rep(one2n - 0.5, p)
        xr <- xl + width
        yt <- yb + 1

        if (!is.null(labx <- colnames(x)))
            text(xx[-length(xx)] + width/2, max(yt), labx, adj = c(0.5, -0.5), xpd = TRUE)

        if (continuous) {
            nux <- if (is.logical(continuous)) 10 else continuous
            sq <- seq(min(x), max(x), length.out = nux + 1)
            x <- .bincode(x, sq, FALSE, TRUE)
            lgd <- paste0("[", sq[-length(sq)], "-", sq[-1], ")")
        } else {
            if (is.raw(x)) x <- toupper(as.character(x)) # for DNAbin objects
            nux <- length(ux <- sort(unique.default(x)))
            x <- match(x, ux)
            lgd <- as.character(ux)
        }
        co <- funcol(nux)
        conames <- names(co)
        if (!is.null(conames)) co <- co[lgd]
        rect(xl, yb, xr, yt, col = co[x], xpd = TRUE, ...)
        legend <- match.arg(legend, c("below", "side", "none"))
        if (legend != "none") {
            if (legend == "below")
                legend((x0 + x1)/2, -yinch(0.1), lgd, pch = 22, pt.bg = co,
                       pt.cex = 2, bty = "n", xjust = 0.5, yjust = 0.5,
                       horiz = TRUE, xpd = TRUE)
            else legend(x1, n, lgd, pch = 22, pt.bg = co,
                        pt.cex = 2, bty = "n", yjust = 1, xpd = TRUE)
        }
    })
}
