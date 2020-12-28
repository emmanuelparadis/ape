## makeLabel.R (2019-10-14)

##   Label Management

## Copyright 2010-2019 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

makeLabel <- function(x, ...) UseMethod("makeLabel")

makeLabel.character <- function(x, len = 99, space = "_",
          make.unique = TRUE, illegal = "():;,[]", quote = FALSE, ...)
{
    x <- gsub("[[:space:]]", space, x)
    if (illegal != "") {
        illegal <- unlist(strsplit(illegal, NULL))
        for (i in illegal) x <- gsub(i, "", x, fixed = TRUE)
    }
    if (quote) len <- len - 2
    nc <- nchar(x) > len
    if (any(nc)) x[nc] <- substr(x[nc], 1, len)
    tab <- table(x)
    if (all(tab == 1)) make.unique <- FALSE
    if (make.unique) {
        dup <- tab[which(tab > 1)]
        nms <- names(dup)
        for (i in 1:length(dup)) {
            j <- which(x == nms[i])
            end <- nchar(x[j][1])
            ## w: number of characters to be added as suffix
            w <- floor(log10(dup[i])) + 1
            suffix <- formatC(1:dup[i], width = w, flag = "0")
            if (end + w > len) {
                start <- end - w + 1
                substr(x[j], start, end) <- suffix
            } else x[j] <- paste(x[j], suffix, sep = "")
        }
    }
    if (quote) x <- paste('"', x, '"', sep = "")
    x
}

makeLabel.phylo <- function(x, tips = TRUE, nodes = TRUE, ...)
{
    if (tips)
        x$tip.label <- makeLabel.character(x$tip.label, ...)
    if (!is.null(x$node.label) && nodes)
        x$node.label <- makeLabel.character(x$node.label, ...)
    x
}

makeLabel.multiPhylo <- function(x, tips = TRUE, nodes = TRUE, ...)
{
    y <- attr(x, "TipLabel")
    if (is.null(y)) {
        for (i in 1:length(x))
            x[[i]] <- makeLabel.phylo(x[[i]], tips = tips, nodes = nodes, ...)
    } else {
        attr(x, "TipLabel") <- makeLabel.character(y, ...)
    }
    x
}

makeLabel.DNAbin <- function(x, ...)
{
    if (is.list(x))
        names(x) <- makeLabel.character(names(x), ...)
    else rownames(x) <- makeLabel.character(rownames(x), ...)
    x
}

mixedFontLabel <-
    function(..., sep = " ", italic = NULL, bold = NULL, parenthesis = NULL,
             always.upright = c("sp.", "spp.", "ssp."))
{
    x <- list(...)
    n <- length(x)

    if (!is.null(italic)) {
        for (i in italic) {
            y <- x[[i]]
            s <- ! y %in% always.upright
            y[s] <- paste("italic(\"", y[s], "\")", sep = "")
            if (any(!s)) y[!s] <- paste("plain(\"", y[!s], "\")", sep = "")
            x[[i]] <- y
        }
    }

    if (!is.null(bold)) {
        for (i in bold) {
            y <- x[[i]]
            s <- logical(length(y))
            s[grep("^italic", y)] <- TRUE
            y[s] <- sub("^italic", "bolditalic", y[s])
            y[!s] <- paste("bold(\"", y[!s], "\")", sep = "")
            x[[i]] <- y
        }
    }

    k <- which(! 1:n %in% c(italic, bold)) # those in upright
    if (length(k))
        for (i in k)
            x[[i]] <- paste("plain(\"", x[[i]], "\")", sep = "")

    if (!is.null(parenthesis))
        for (i in parenthesis)
            x[[i]] <- paste("(", x[[i]], ")", sep = "")

    res <- x[[1L]]
    if (n > 1) {
        sep <- rep(sep, length.out = n - 1L)
        for (i in 2:n)
            res <- paste(res, "*\"", sep[i - 1L], "\"*", x[[i]], sep = "")
    }
    parse(text = res)
}

.getSeparatorTaxaLabels <- function(x)
{
    if (length(grep("_", x))) "_" else " "
}

label2table <- function(x, sep = NULL, as.is = FALSE)
{
    n <- length(x)
    if (is.null(sep)) sep <- .getSeparatorTaxaLabels(x)
    x <- strsplit(x, sep)
    maxlen <- max(lengths(x))
    x <- unlist(lapply(x, "[", 1:maxlen))
    x <- matrix(x, n, maxlen, byrow = TRUE)
    x <- as.data.frame(x, as.is = as.is)
    baselevels <- c("genus", "species", "subspecies")
    nmx <- if (maxlen <= 3) baselevels[1:maxlen]
           else c(baselevels, paste0("type", 1:(maxlen - 3)))
    names(x) <- nmx
    x
}

stripLabel <- function(x, species = FALSE, subsp = TRUE, sep = NULL)
{
    if (is.null(sep)) sep <- .getSeparatorTaxaLabels(x)
    n <- 0
    if (species) n <- 1 else if (subsp) n <- 2
    if (!n) return(x)
    x <- strsplit(x, sep)
    x <- lapply(x, "[", 1:n)
    sapply(x, paste, collapse = sep)
}

abbreviateGenus <- function(x, genus = TRUE, species = FALSE, sep = NULL)
{
    if (is.null(sep)) sep <- .getSeparatorTaxaLabels(x)
    if (genus) x <- sub(paste0("[[:lower:]]{1,}", sep), paste0(".", sep), x)
    if (!species) return(x)
    x <- strsplit(x, sep)
    k <- which(lengths(x, use.names = FALSE) > 1)
    for (i in k)
        x[[i]][2] <- paste0(substr(x[[i]][2], 1, 1), ".")
    sapply(x, paste, collapse = sep)
}

updateLabel <- function(x, old, new, ...) UseMethod("updateLabel")

updateLabel.character <- function(x, old, new, exact = TRUE, ...)
{
    if (length(old) != length(new))
        stop("'old' and 'new' not of the same length")
    if (exact) {
        for (i in seq_along(old))
            x[x == old[i]] <- new[i]
    } else {
        for (i in seq_along(old))
            x[grep(old[i], x)] <- new[i]
    }
    x
}

updateLabel.DNAbin <- function(x, old, new, exact = TRUE, ...)
{
    labs <- labels(x)
    labs <- updateLabel.character(labs, old, new, exact, ...)
    if (is.list(x)) names(x) <- labs else rownames(x) <- labs
    x
}

updateLabel.AAbin <- function(x, old, new, exact = TRUE, ...)
    updateLabel.DNAbin(x, old, new, exact, ...)

updateLabel.phylo <- function(x, old, new, exact = TRUE, nodes = FALSE, ...)
{
    x$tip.label <- updateLabel.character(x$tip.label, old, new, exact, ...)
    if (nodes)
        x$node.label <- updateLabel.character(x$node.label, old, new, exact, ...)
    x
}

updateLabel.evonet <- function(x, old, new, exact = TRUE, nodes = FALSE, ...)
    updateLabel.phylo(x, old, new, exact, nodes, ...)


updateLabel.data.frame <- function(x, old, new, exact = TRUE, ...)
{
    row.names(x) <- updateLabel.character(row.names(x), old, new, exact, ...)
    x
}

updateLabel.matrix <- function(x, old, new, exact = TRUE, ...)
{
    rownames(x) <- updateLabel.character(rownames(x), old, new, exact, ...)
    x
}
