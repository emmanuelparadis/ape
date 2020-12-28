## nodelabels.R (2020-12-06)

##   Labelling Trees

## Copyright 2004-2020 Emmanuel Paradis, 2006 Ben Bolker, and 2006 Jim Lemon

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

## from JL:
## floating.pie() from plotrix with two changes:
## (1) aspect ratio fixed, so pies will appear circular
##     (`radius' is the radius in user coordinates along the x axis);
## (2) zero values allowed (but not negative).

floating.pie.asp <- function(xpos, ypos, x, edges = 200, radius = 1,
                             col = NULL, startpos = 0, ...)
{
    u <- par("usr")
    user.asp <- diff(u[3:4])/diff(u[1:2])
    p <- par("pin")
    inches.asp <- p[2]/p[1]
    asp <- user.asp/inches.asp
    if (!is.numeric(x) || any(is.na(x) | x < 0))
        stop("floating.pie: x values must be non-negative")
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    col <- if (is.null(col)) rainbow(nx) else rep_len(col, nx)
    ## next a fix from Klaus to avoid a "3-o'clock" segment on pies with
    ## only one proportion equal to 1:
    if (length(i <- which(dx == 1))) {
        symbols(xpos, ypos, circles = radius, inches = FALSE, add = TRUE,
                fg = par("fg"), bg = col[i]) # suggested by Liam
    } else {
        bc <- 2 * pi * (x[1:nx] + dx/2) + startpos
        for (i in seq_len(nx)) {
            n <- max(2, floor(edges * dx[i]))
            t2p <- 2 * pi * seq(x[i], x[i + 1], length = n) + startpos
            xc <- c(cos(t2p) * radius + xpos, xpos)
            yc <- c(sin(t2p) * radius*asp + ypos, ypos)
            polygon(xc, yc, col = col[i], ...)
        }
    }
}

BOTHlabels <- function(text, sel, XX, YY, adj, frame, pch, thermo,
                       pie, piecol, col, bg, horiz, width, height, ...)
{
    if (missing(text)) text <- NULL
    if (length(adj) == 1) adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie))
        text <- as.character(sel)
    frame <- match.arg(frame, c("rect", "circle", "none"))
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) args$cex else par("cex")
    if (frame != "none" && !is.null(text)) {
        if (frame == "rect") {
            width <- strwidth(text, units = "inches", cex = CEX)
            height <- strheight(text, units = "inches", cex = CEX)
            if ("srt" %in% names(args)) {
                args$srt <- args$srt %% 360 # just in case srt >= 360
                if (args$srt == 90 || args$srt == 270) {
                    tmp <- width
                    width <- height
                    height <- tmp
                } else if (args$srt != 0)
                  warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
            }
            width <- xinch(width)
            height <- yinch(height)
            xl <- XX - width*adj[1] - xinch(0.03)
            xr <- xl + width + xinch(0.03)
            yb <- YY - height*adj[2] - yinch(0.02)
            yt <- yb + height + yinch(0.05)
            rect(xl, yb, xr, yt, col = bg)
        }
        if (frame == "circle") {
            radii <- 0.8*apply(cbind(strheight(text, units = "inches", cex = CEX),
                                     strwidth(text, units = "inches", cex = CEX)), 1, max)
            symbols(XX, YY, circles = radii, inches = max(radii), add = TRUE, bg = bg)
        }
    }
    if (!is.null(thermo)) {
        parusr <- par("usr")

        if (is.null(width)) {
            width <- CEX * (parusr[2] - parusr[1])
            width <- if (horiz) width/15 else width/40
        }

        if (is.null(height)) {
            height <- CEX * (parusr[4] - parusr[3])
            height <- if (horiz) height/40 else height/15
        }

        if (is.vector(thermo) || ncol(thermo) == 1) thermo <- cbind(thermo, 1 - thermo)
        thermo <- if (horiz) width * thermo else height * thermo
        if (is.null(piecol)) piecol <- rainbow(ncol(thermo))

        xl <- XX - width/2 + adj[1] - 0.5 # added 'adj' from Janet Young (2009-09-30)
        xr <- xl + width
        yb <- YY - height/2 + adj[2] - 0.5
        yt <- yb + height

        if (horiz) {
            ## draw the first rectangle:
            rect(xl, yb, xl + thermo[, 1], yt, border = NA, col = piecol[1])
            for (i in 2:ncol(thermo))
                rect(xl + rowSums(thermo[, 1:(i - 1), drop = FALSE]), yb,
                     xl + rowSums(thermo[, 1:i]), yt, border = NA, col = piecol[i])
        } else {
            ## draw the first rectangle:
            rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
            for (i in 2:ncol(thermo))
                rect(xl, yb + rowSums(thermo[, 1:(i - 1), drop = FALSE]),
                     xr, yb + rowSums(thermo[, 1:i]),
                     border = NA, col = piecol[i])
        }

        ## check for NA's before drawing the borders
        s <- apply(thermo, 1, function(xx) any(is.na(xx)))
        xl[s] <-  xr[s] <- NA
        rect(xl, yb, xr, yt, border = "black")

        if (!horiz) {
            segments(xl, YY, xl - width/5, YY)
            segments(xr, YY, xr + width/5, YY)
        }
    }
    ## from BB:
    if (!is.null(pie)) {
        if (is.data.frame(pie)) pie <- as.matrix(pie)
        if (is.vector(pie) || ncol(pie) == 1) pie <- cbind(pie, 1 - pie)
        xrad <- CEX * diff(par("usr")[1:2]) / 50
        xrad <- rep(xrad, length(sel))
        XX <- XX + adj[1] - 0.5
        YY <- YY + adj[2] - 0.5
        for (i in seq_along(sel)) {
            if (any(is.na(pie[i, ]))) next
            floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], col = piecol)
        }
    }
    if (!is.null(text)) text(XX, YY, text, adj = adj, col = col, ...)
    if (!is.null(pch)) points(XX + adj[1] - 0.5, YY + adj[2] - 0.5,
                              pch = pch, col = col, bg = bg, ...)
}

nodelabels <-
    function(text, node, adj = c(0.5, 0.5), frame = "rect",
             pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
             col = "black", bg = "lightblue", horiz = FALSE,
             width = NULL, height = NULL, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(node)) node <- (lastPP$Ntip + 1):length(lastPP$xx)
    XX <- lastPP$xx[node]
    YY <- lastPP$yy[node]
    BOTHlabels(text, node, XX, YY, adj, frame, pch, thermo,
               pie, piecol, col, bg, horiz, width, height, ...)
}

tiplabels <-
    function(text, tip, adj = c(0.5, 0.5), frame = "rect",
             pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
             col = "black", bg = "yellow", horiz = FALSE,
             width = NULL, height = NULL, offset = 0, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(tip)) tip <- 1:lastPP$Ntip
    XX <- lastPP$xx[tip]
    YY <- lastPP$yy[tip]
    if (offset != 0) {
        if (lastPP$type %in% c("phylogram", "cladogram")) {
            switch(lastPP$direction,
                   "rightwards" = {XX <- XX + offset},
                   "leftwards" = {XX <- XX - offset},
                   "upwards" = {YY <- YY + offset},
                   "downwards" = {YY <- YY - offset})
        } else {
            if (lastPP$type %in% c("fan", "radial")) {
                tmp <- rect2polar(XX, YY)
                if (lastPP$align.tip.label) tmp$r[] <- max(tmp$r)
                tmp <- polar2rect(tmp$r + offset, tmp$angle)
                XX <- tmp$x
                YY <- tmp$y
            } else {
                if (lastPP$type == "unrooted")
                    warning("argument 'offset' ignored with unrooted trees")
            }
        }
    }
    BOTHlabels(text, tip, XX, YY, adj, frame, pch, thermo,
               pie, piecol, col, bg, horiz, width, height, ...)
}

edgelabels <-
    function(text, edge, adj = c(0.5, 0.5), frame = "rect",
             pch = NULL, thermo = NULL, pie = NULL, piecol = NULL,
             col = "black", bg = "lightgreen", horiz = FALSE,
             width = NULL, height = NULL, date = NULL, ...)
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(edge)) {
        sel <- 1:dim(lastPP$edge)[1]
        subedge <- lastPP$edge
    } else {
        sel <- edge
        subedge <- lastPP$edge[sel, , drop = FALSE]
    }

    xx <- lastPP$xx
    yy <- lastPP$yy

    if (lastPP$type == "phylogram") {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            XX <- (xx[subedge[, 1]] + xx[subedge[, 2]]) / 2
            YY <- yy[subedge[, 2]]
        } else {
            XX <- xx[subedge[, 2]]
            YY <- (yy[subedge[, 1]] + yy[subedge[, 2]]) / 2
        }
    } else {
        if (lastPP$type == "fan") { # fix by Klaus Schliep (2015-07-31)
            r <- sqrt(xx^2 + yy^2)
            tmp <- (r[subedge[, 2]] + r[subedge[, 1]]) / (r[subedge[, 2]] * 2)
            XX <- xx[subedge[, 2]] * tmp
            YY <- yy[subedge[, 2]] * tmp
        } else {
            XX <- (xx[subedge[, 1]] + xx[subedge[, 2]]) / 2
            YY <- (yy[subedge[, 1]] + yy[subedge[, 2]]) / 2
        }
    }

    ## suggestion by Rob Lanfear:
    if (!is.null(date)) XX[] <- max(lastPP$xx) - date

    BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo,
               pie, piecol, col, bg, horiz, width, height, ...)
}

edges <- function(nodes0, nodes1, arrows = 0, type = "classical", ...)
{
    type <- match.arg(type, c("classical", "triangle", "harpoon"))
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    ## we do the recycling if necessary:
    if (length(nodes0) != length(nodes1)) {
        tmp <- cbind(nodes0, nodes1)
        nodes0 <- tmp[, 1]
        nodes1 <- tmp[, 2]
    }
    x0 <- lastPP$xx[nodes0]
    y0 <- lastPP$yy[nodes0]
    x1 <- lastPP$xx[nodes1]
    y1 <- lastPP$yy[nodes1]
    if (arrows)
        if (type == "classical")
            arrows(x0, y0, x1, y1, code = arrows, ...)
        else
            fancyarrows(x0, y0, x1, y1, code = arrows, type = type, ...)
    else
        segments(x0, y0, x1, y1, ...)
}

fancyarrows <-
    function(x0, y0, x1, y1, length = 0.25, angle = 30, code = 2,
             col = par("fg"), lty = par("lty"), lwd = par("lwd"),
             type = "triangle", ...)
{
    foo <- function(x0, y0, x1, y1) {
        ## important to correct with these parameters cause
        ## the coordinate system will likely not be Cartesian
        pin <- par("pin")
        usr <- par("usr")
        A1 <- pin[1]/diff(usr[1:2])
        A2 <- pin[2]/diff(usr[3:4])
        x0 <- x0 * A1
        y0 <- y0 * A2
        x1 <- x1 * A1
        y1 <- y1 * A2
        atan2(y1 - y0, x1 - x0)
    }
    arrow.triangle <- function(x, y) {
        beta <- alpha - angle/2
        xa <- xinch(length * cos(beta)) + x
        ya <- yinch(length * sin(beta)) + y
        beta <- beta + angle
        xb <- xinch(length * cos(beta)) + x
        yb <- yinch(length * sin(beta)) + y
        n <- length(x)
        col <- rep(col, length.out = n)
        for (i in 1:n)
            polygon(c(x[i], xa[i], xb[i]), c(y[i], ya[i], yb[i]),
                    col = col[i], border = col[i])
        list((xa + xb)/2, (ya + yb)/2)
    }
    arrow.harpoon <- function(x, y) {
        beta <- alpha - angle/2
        xa <- xinch(length * cos(beta)) + x
        ya <- yinch(length * sin(beta)) + y
        beta <- alpha + angle/2
        xb <- xinch(length * cos(beta)) + x
        yb <- yinch(length * sin(beta)) + y
        xc <- x/2 + (xa + xb)/4
        yc <- y/2 + (ya + yb)/4
        n <- length(x)
        col <- rep(col, length.out = n)
        for (i in 1:n)
            polygon(c(x[i], xa[i], xc[i], xb[i]),
                    c(y[i], ya[i], yc[i], yb[i]),
                    col = col[i], border = col[i])
        list(xc, yc)
    }

    type <- match.arg(type, c("triangle", "harpoon"))
    angle <- pi*angle/180 # degree -> radian
    alpha <- foo(x0, y0, x1, y1) # angle of segment with x-axis
    ## alpha is in [-pi, pi]

    FUN <- if (type == "triangle") arrow.triangle else arrow.harpoon
    XY0 <- if (code == 1 || code == 3) FUN(x0, y0) else list(x0, y0)
    if (code >= 2) {
        alpha <- (alpha + pi) %% (2 * pi)
        XY1 <- FUN(x1, y1)
    } else XY1 <- list(x1, y1)
    segments(XY0[[1]], XY0[[2]], XY1[[1]], XY1[[2]], col = col, lty = lty, lwd = lwd, ...)
}
