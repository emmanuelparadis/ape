## read.tree.R (2021-05-04)

##   Read Tree Files in Parenthetic Format

## Copyright 2002-2021 Emmanuel Paradis, Daniel Lawson and Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

read.tree <- function(file = "", text = NULL, tree.names = NULL, skip = 0,
    comment.char = "", keep.multi = FALSE, ...)
{
    if (!is.null(text)) {
        if (!is.character(text))
            stop("argument 'text' must be of mode character")
        tree <- text
    } else {
        tree <- scan(file = file, what = "", sep = "\n", quiet = TRUE,
                     skip = skip, comment.char = comment.char, ...)
    }

    ## Suggestion from Eric Durand and Nicolas Bortolussi (added 2005-08-17):
    if (identical(tree, character(0))) {
        warning("empty character string.")
        return(NULL)
    }

#    tree <- gsub("'", "", tree) # fixed by GuangchuangYu (2021-01-04)

    single_quotes <- function(x) {
        single.quote <- as.raw(39)
        x <- charToRaw(x)
        z <- which(x == single.quote)
        if (length(z) %% 2)
            stop("wrong number of single quotes around labels")
        l <- length(z) / 2
        opening <- z[c(TRUE, FALSE)]
        closing <- z[c(FALSE, TRUE)]
        from <- c(1, closing + 1L)
        to <- c(opening - 1L, length(x))
        i <- mapply(":", from = from, to = to, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        keep <- lapply(i, function(i) x[i])
        tmp_label <- paste0("@_", 1:l, "_@")
        tmpLabsRaw <- lapply(tmp_label, charToRaw)
        n <- 2 * l + 1L
        res <- vector("list", n)
        res[seq(1, n, 2)] <- keep
        res[seq(2, n - 1, 2)] <- tmpLabsRaw
        res <- rawToChar(unlist(res))
        i <- mapply(":", from = opening, to = closing, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        orig_label <- lapply(i, function(i) x[i])
        orig_label <- sapply(orig_label, rawToChar)
        names(orig_label) <- tmp_label
        list(res, orig_label)
    }

    y <- unlist(gregexpr(";", tree))

    ## if one tree per line much faster
    if (identical(y, nchar(tree))) { # check if always one tree per line
        Ntree <- length(y)
        STRING <- character(Ntree)
        for (i in 1:Ntree) {
            STRING[i] <- gsub("\\[[^]]*\\]", "", tree[i]) # delete comments (fix 2015-01-12)
        }
    } else {
        tree <- unlist(strsplit(tree, NULL))
        y <- which(tree == ";")
        Ntree <- length(y)
        x <- c(1, y[-Ntree] + 1)
        ## Suggestion from Olivier Francois (added 2006-07-15):
        if (is.na(y[1])) return(NULL)
        STRING <- character(Ntree)
        for (i in 1:Ntree) {
            tmp <- paste0(tree[x[i]:y[i]], collapse = "")
            STRING[i] <- gsub("\\[[^]]*\\]", "", tmp) # delete comments (fix 2015-01-12)
        }
    }

    ## remove possible leading and trailing underscores
    STRING <- gsub("^_+", "", STRING)
    STRING <- gsub("_+$", "", STRING)

    ## replace labels with single quotes (moved from above)
    z <- grepl("'", STRING)
    if (any(z)) {
#        Ntree <- length(tree)
        tmp_label <- vector("list", Ntree)
        for (i in 1:Ntree) {
            if (z[i]) {
                TMP <- single_quotes(STRING[i])
                STRING[i] <- TMP[[1]]
                tmp_label[[i]] <- TMP[[2]]
            }
        }
    }

    tree <- gsub("[ \t]", "", tree) # moved from above: now spaces within (single) quoted labels are not deleted

    getTreeName <- function(x) {
        res <- rep("", length(x))
        i <- regexpr("\\(", x)
        s <- i > 1
        if (any(s)) res[s] <- substr(x[s], 1, i[s] - 1)
        res
    }

    tmpnames <- getTreeName(STRING)
    if (is.null(tree.names) && any(nzchar(tmpnames))) tree.names <- tmpnames

    colon <- grep(":", STRING)
    if (!length(colon)) {
        obj <- lapply(STRING, .cladoBuild)
    } else if (length(colon) == Ntree) {
        obj <- lapply(STRING, .treeBuild)
    } else {
        obj <- vector("list", Ntree)
        obj[colon] <- lapply(STRING[colon], .treeBuild)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        obj[nocolon] <- lapply(STRING[nocolon], .cladoBuild)
    }

    for (i in 1:Ntree) {
        if (z[i]) {
            #browser()
            tmp_lab <- tmp_label[[i]]
            tip.label <- obj[[i]]$tip.label
            node.label <- obj[[i]]$node.label
            ind <- match(tip.label, names(tmp_lab))
            ind2 <- which(!is.na(ind))
            if (length(ind2)) {
                tip.label[ind2] <- tmp_lab[ind[ind2]]
                #tmp_lab <- tmp_lab[-ind[ind2]]
            }

            ind <- match(node.label, names(tmp_lab))
            ind2 <- which(!is.na(ind))
            if (length(ind2)) {
                node.label[ind2] <- tmp_lab[ind[ind2]]
                #tmp_lab <- tmp_lab[-ind[ind2]]
            }
            #if (length(tmp_lab)) {
            #    for (j in 1:length(tmp_lab)) {
            #        node.label <- gsub(names(tmp_lab)[j], tmp_lab[j], node.label)
            #        tip.label <- gsub(names(tmp_lab)[j], tmp_lab[j], tip.label)
            #    }
            #}
            obj[[i]]$tip.label <- tip.label
            obj[[i]]$node.label <- node.label
        }
    }
    if (Ntree == 1 && !keep.multi) obj <- obj[[1]] else {
        if (!is.null(tree.names)) names(obj) <- tree.names
        class(obj) <- "multiPhylo"
    }
    obj
}
