## cophyloplot.R (2014-04-07)

##   Plots two phylogenetic trees face to
##   face with the links between the tips

## Copyright 2008-2010 Damien de Vienne

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

cophyloplot <-
    function(x, y, assoc = NULL, use.edge.length = FALSE, space = 0,
             length.line = 1, gap = 2, type = "phylogram", rotate = FALSE,
             col = par("fg"), lwd = par("lwd"), lty = par("lty"),
             show.tip.label = TRUE, font = 3, ...)
{
    if (is.null(assoc)) {
        assoc <- matrix(ncol = 2)
        print("No association matrix specified. Links will be omitted.")
    }
    if (rotate == TRUE) {
        cat("\n    Click on a node to rotate (right click to exit)\n\n")
        repeat {
            res <- plotCophylo2(x, y, assoc = assoc, use.edge.length = use.edge.length,
                space = space, length.line = length.line, gap = gap,
                type = type, return = TRUE, col = col, lwd=lwd, lty=lty, show.tip.label = show.tip.label,
                font = font)
            click <- identify(res$c[, 1], res$c[, 2], n = 1)
            if (click < length(res$a[, 1]) + 1) {
                if (click > res$N.tip.x)
                  x <- rotate(x, click)
            } else if (click < length(res$c[, 1]) + 1) {
                if (click > length(res$a[, 1]) + res$N.tip.y)
                  y <- rotate(y, click - length(res$a[, 1]))
            }
        }
        on.exit(cat("done\n"))
    }
    else plotCophylo2(x, y, assoc = assoc, use.edge.length = use.edge.length,
        space = space, length.line = length.line, gap = gap,
        type = type, return = FALSE, col = col, lwd=lwd, lty=lty, show.tip.label = show.tip.label, font = font)
}

plotCophylo2 <-
    function(x, y, assoc = assoc, use.edge.length = use.edge.length,
             space = space, length.line = length.line, gap = gap,
             type = type, return = return, col = col, lwd=lwd, lty=lty,
             show.tip.label = show.tip.label,
             font = font, ...)
{
    res <- list()
###choice of the minimum space between the trees
    left <- max(nchar(x$tip.label, type = "width")) + length.line
    right <- max(nchar(y$tip.label, type = "width")) + length.line
    space.min <- left + right + gap * 2
    if ((space <= 0) || (space < space.min)) space <- space.min
    N.tip.x <- Ntip(x)
    N.tip.y <- Ntip(y)
    res$N.tip.x <- N.tip.x
    res$N.tip.y <- N.tip.y
    a <- plotPhyloCoor(x, use.edge.length = use.edge.length, type = type)
    res$a <- a
    b <- plotPhyloCoor(y, use.edge.length = use.edge.length,
                       direction = "leftwards", type = type)
###for the two trees to have the extreme leaves at the same ordinate.
    a[, 2] <- a[, 2] - min(a[, 2])
    b[, 2] <- b[, 2] - min(b[, 2])
    res$b <- b
    b2 <- b
    b2[, 1] <- b[1:nrow(b), 1] * (max(a[, 1])/max(b[, 1])) +
        space + max(a[, 1])
    b2[, 2] <- b[1:nrow(b), 2] * (max(a[, 2])/max(b[, 2]))
    res$b2 <- b2
    c <- matrix(ncol = 2, nrow = nrow(a) + nrow(b))
    c[1:nrow(a), ] <- a[1:nrow(a), ]
    c[nrow(a) + 1:nrow(b), 1] <- b2[, 1]
    c[nrow(a) + 1:nrow(b), 2] <- b2[, 2]
    res$c <- c
    plot(c, type = "n", xlim = NULL, ylim = NULL, log = "", main = NULL,
        sub = NULL, xlab = NULL, ylab = NULL, ann = FALSE, axes = FALSE,
        frame.plot = FALSE)
 ###segments for cladograms
   if (type == "cladogram") {
        for (i in 1:(nrow(a) - 1)) segments(a[x$edge[i, 1], 1],
            a[x$edge[i, 1], 2], a[x$edge[i, 2], 1], a[x$edge[i,
                2], 2], col="red")
        for (i in 1:(nrow(b) - 1))
            segments(b2[y$edge[i, 1], 1], b2[y$edge[i, 1], 2],
                     b2[y$edge[i, 2], 1], b2[y$edge[i, 2], 2])
    }
###segments for phylograms
    if (type == "phylogram") {
        for (i in (N.tip.x + 1):nrow(a)) {
            l <- length(x$edge[x$edge[, 1] == i, ][, 1])
            for (j in 1:l) {
                segments(a[x$edge[x$edge[, 1] == i, ][1, 1],
                  1], a[x$edge[x$edge[, 1] == i, 2], 2][1], a[x$edge[x$edge[,
                  1] == i, ][1, 1], 1], a[x$edge[x$edge[, 1] ==
                  i, 2], 2][j])
                segments(a[x$edge[x$edge[, 1] == i, ][1, 1],
                  1], a[x$edge[x$edge[, 1] == i, 2], 2][j], a[x$edge[x$edge[,
                  1] == i, 2], 1][j], a[x$edge[x$edge[, 1] ==
                  i, 2], 2][j])
            }
        }
        for (i in (N.tip.y + 1):nrow(b)) {
            l <- length(y$edge[y$edge[, 1] == i, ][, 1])
            for (j in 1:l) {
                segments(b2[y$edge[y$edge[, 1] == i, ][1, 1],
                  1], b2[y$edge[y$edge[, 1] == i, 2], 2][1],
                  b2[y$edge[y$edge[, 1] == i, ][1, 1], 1], b2[y$edge[y$edge[,
                    1] == i, 2], 2][j])
                segments(b2[y$edge[y$edge[, 1] == i, ][1, 1],
                  1], b2[y$edge[y$edge[, 1] == i, 2], 2][j],
                  b2[y$edge[y$edge[, 1] == i, 2], 1][j], b2[y$edge[y$edge[,
                    1] == i, 2], 2][j])
            }
        }
    }
    if (show.tip.label) {
        text(a[1:N.tip.x, ], cex = 0, font = font, pos = 4,
             labels = x$tip.label)
        text(b2[1:N.tip.y, ], cex = 1, font = font, pos = 2,
             labels = y$tip.label)
    }
###links between associated taxa. Takes into account the size of the character strings of the taxa names.
    lsa <- 1:N.tip.x
    lsb <- 1:N.tip.y
    decx <- array(nrow(assoc))
    decy <- array(nrow(assoc))

    #colors
    if (length(col)==1) colors<-c(rep(col, nrow(assoc)))
    else if (length(col)>=nrow(assoc)) colors<-col
    else  colors<-c(rep(col, as.integer(nrow(assoc)/length(col))+1))

    #lwd
    if (length(lwd)==1) lwidths<-c(rep(lwd, nrow(assoc)))
    else if (length(lwd)>=nrow(assoc)) lwidths<-lwd
    else  lwidths<-c(rep(lwd, as.integer(nrow(assoc)/length(lwd))+1))

    #lty
    if (length(lty) == 1) ltype <- c(rep(lty, nrow(assoc)))
    else if (length(lty) >= nrow(assoc)) ltype <- lty
    else  ltype <- c(rep(lty, as.integer(nrow(assoc)/length(lty))+1))

    for (i in 1:nrow(assoc)) {
        if (show.tip.label) {
            decx[i] <- strwidth(x$tip.label[lsa[x$tip.label == assoc[i, 1]]])
            decy[i] <- strwidth(y$tip.label[lsb[y$tip.label == assoc[i, 2]]])
        } else {
            decx[i] <- decy[i] <- 0
        }

        if (length.line) { # added by EP (2014-04-07)
            segments(a[lsa[x$tip.label == assoc[i, 1]], 1] + decx[i] + gap,
                     a[lsa[x$tip.label == assoc[i, 1]], 2],
                     a[lsa[x$tip.label == assoc[i, 1]], 1] + gap + left,
                     a[lsa[x$tip.label ==  assoc[i, 1]], 2],
                     col = colors[i], lwd = lwidths[i], lty = ltype[i])

            segments(b2[lsb[y$tip.label == assoc[i, 2]], 1] - (decy[i] + gap),
                     b2[lsb[y$tip.label == assoc[i, 2]], 2],
                     b2[lsb[y$tip.label == assoc[i, 2]], 1] - (gap + right),
                     b2[lsb[y$tip.label ==  assoc[i, 2]], 2],
                     col = colors[i], lwd = lwidths[i], lty = ltype[i])
        }

        segments(a[lsa[x$tip.label == assoc[i, 1]], 1] + gap + left,
                 a[lsa[x$tip.label == assoc[i, 1]], 2],
                 b2[lsb[y$tip.label == assoc[i, 2]], 1] - (gap + right),
                 b2[lsb[y$tip.label == assoc[i, 2]], 2],
                 col = colors[i], lwd = lwidths[i], lty = ltype[i])
    }
    if (return == TRUE)  return(res)
}
