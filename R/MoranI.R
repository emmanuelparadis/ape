## MoranI.R (2008-01-14)

##   Moran's I Autocorrelation Index

## Copyright 2004 Julien Dutheil, 2007-2008 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

## code cleaned-up by EP (Dec. 2007)

Moran.I <- function(x, weight, scaled = FALSE, na.rm = FALSE,
                    alternative = "two.sided")
{
    if(dim(weight)[1] != dim(weight)[2])
        stop("'weight' must be a square matrix")
    n <- length(x)
    if(dim(weight)[1] != n)
        stop("'weight' must have as many rows as observations in 'x'")
    ## Expected mean:
    ei <- -1/(n - 1)

    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        } else {
            warning("'x' has missing values: maybe you wanted to set na.rm = TRUE?")
            return(list(observed = NA, expected = ei, sd = NA, p.value = NA))
        }
    }

    ## normalizing the weights:
    ## Note that we normalize after possibly removing the
    ## missing data.
    ROWSUM <- rowSums(weight)
    ## the following is useful if an observation has no "neighbour":
    ROWSUM[ROWSUM == 0] <- 1
    weight <- weight/ROWSUM # ROWSUM is properly recycled

    s <- sum(weight)
    m <- mean(x)
    y <- x - m # centre the x's
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    obs <- (n/s) * (cv/v)
    ## Scaling:
    if (scaled) {
        i.max <- (n/s) * (sd(rowSums(weight) * y)/sqrt(v/(n - 1)))
        obs <- obs/i.max
    }
    ## Expected sd:
    S1 <- 0.5 * sum((weight + t(weight))^2)
    S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
    ## the above is the same than:
    ##S2 <- 0
    ##for (i in 1:n)
    ##    S2 <- S2 + (sum(weight[i, ]) + sum(weight[, i]))^2

    s.sq <- s^2
    k <- (sum(y^4)/n) / (v/n)^2
    sdi <- sqrt((n*((n^2 - 3*n + 3)*S1 - n*S2 + 3*s.sq) -
                 k*(n*(n - 1)*S1 - 2*n*S2 + 6*s.sq))/
                ((n - 1)*(n - 2)*(n - 3)*s.sq) - 1/((n - 1)^2))

    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
    pv <- pnorm(obs, mean = ei, sd = sdi)
    if (alternative == "two.sided")
        pv <- if (obs <= ei) 2*pv else 2*(1 - pv)
    if (alternative == "greater") pv <- 1 - pv
    list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}

weight.taxo <- function(x)
{
    d <- outer(x, x, "==")
    diag(d) <- 0 # implicitly converts 'd' into numeric
    d
}

weight.taxo2 <- function(x, y)
{
    d <- outer(x, x, "==") & outer(y, y, "!=")
    diag(d) <- 0
    d
}

correlogram.formula <- function(formula, data = NULL, use = "all.obs")
{
    err <- 'formula must be of the form "y1+...+yn ~ x1/x2/../xn"'
    use <- match.arg(use, c("all.obs", "complete.obs", "pairwise.complete.obs"))
    if (formula[[1]] != "~") stop(err)

    lhs <- formula[[2]]
    y.nms <- if (length(lhs) > 1)
        unlist(strsplit(as.character(as.expression(lhs)), " \\+ "))
    else as.character(as.expression(lhs))

    rhs <- formula[[3]]
    gr.nms <- if (length(rhs) > 1)
        rev(unlist(strsplit(as.character(as.expression(rhs)), "/")))
    else as.character(as.expression(rhs))

    if (is.null(data)) {
        ## we 'get' the variables in the .GlobalEnv:
        y <- as.data.frame(sapply(y.nms, get))
        gr <- as.data.frame(sapply(gr.nms, get))
    } else {
        y <- data[y.nms]
        gr <- data[gr.nms]
    }
    if (use == "all.obs") {
        na.fail(y)
        na.fail(gr)
    }
    if (use == "complete.obs") {
        sel <- complete.cases(y, gr)
        y <- y[sel]
        gr <- gr[sel]
    }
    na.rm <- use == "pairwise.complete.obs"

    foo <- function(x, gr, na.rm) {
        res <- data.frame(obs = NA, p.values = NA, labels = colnames(gr))
        for (i in 1:length(gr)) {
            sel <- if (na.rm) !is.na(x) & !is.na(gr[, i]) else TRUE
            xx <- x[sel]
            g <- gr[sel, i]
            w <- if (i > 1) weight.taxo2(g, gr[sel, i - 1]) else weight.taxo(g)
            o <- Moran.I(xx, w, scaled = TRUE)
            res[i, 1] <- o$observed
            res[i, 2] <- o$p.value
        }
        ## We need to specify the two classes; if we specify
        ## only "correlogram", 'res' is coerced as a list
        ## (data frames are of class "data.frame" and mode "list")
        structure(res, class = c("correlogram", "data.frame"))
    }

    if (length(y) == 1) foo(y[[1]], gr, na.rm)
    else structure(lapply(y, foo, gr = gr, na.rm = na.rm),
                   names = y.nms, class = "correlogramList")
}

plot.correlogram <-
    function(x, legend = TRUE, test.level = 0.05,
             col = c("grey", "red"), type = "b", xlab = "",
             ylab = "Moran's I", pch = 21, cex = 2, ...)
{
    BG <- col[(x$p.values < test.level) + 1]
    if (pch > 20 && pch < 26) {
        bg <- col
        col <- CO <- "black"
    } else {
        CO <- BG
        BG <- bg <- NULL
    }
    plot(1:length(x$obs), x$obs, type = type, xaxt = "n", xlab = xlab,
         ylab = ylab, col = CO, bg = BG, pch = pch, cex = cex, ...)
    axis(1, at = 1:length(x$obs), labels = x$labels)
    if (legend)
        legend("top", legend = paste(c("P >=", "P <"), test.level),
               pch = pch, col = col, pt.bg = bg, pt.cex = cex, horiz = TRUE)
}

plot.correlogramList <-
    function(x, lattice = TRUE, legend = TRUE,
             test.level = 0.05, col = c("grey", "red"),
             xlab = "", ylab = "Moran's I",
             type = "b", pch = 21, cex = 2, ...)
{
    n <- length(x)
    obs <- unlist(lapply(x, "[[", "obs"))
    pval <- unlist(lapply(x, "[[", "p.values"))
    gr <- factor(unlist(lapply(x, "[[", "labels")),
                 ordered = TRUE, levels = x[[1]]$labels)
    vars <- gl(n, nlevels(gr), labels = names(x))
    BG <- col[(pval < test.level) + 1]
    if (lattice) {
        ## trellis.par.set(list(plot.symbol=list(pch=19)))
        xyplot(obs ~ gr | vars, xlab = xlab, ylab = ylab,
               panel = function(x, y) {
                   panel.lines(x, y, lty = 2)
                   panel.points(x, y, cex = cex, pch = 19, col = BG)
                   ##lattice::panel.abline(h = 0, lty = 3)
               })
    } else {
        if (pch > 20 && pch < 26) {
            bg <- col
            CO <- rep("black", length(obs))
            col <- "black"
        } else {
            CO <- BG
            BG <- bg <- NULL
        }
        plot(as.numeric(gr), obs, type = "n", xlab = xlab,
             ylab = ylab, xaxt = "n")
        for (i in 1:n) {
            sel <- as.numeric(vars) == i
            lines(as.numeric(gr[sel]), obs[sel], type = type, lty = i,
                  col = CO[sel], bg = BG[sel], pch = pch, cex = cex, ...)
        }
        axis(1, at = 1:length(x[[i]]$obs), labels = x[[i]]$labels)
        if (legend) {
            legend("topright", legend = names(x), lty = 1:n, bty = "n")
            legend("top", legend = paste(c("P >=", "P <"), test.level),
                   pch = pch, col = col, pt.bg = bg, pt.cex = cex, horiz = TRUE)
        }
    }
}
