## read.tree.R (2021-05-26)

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

    ## make a single string
    if (length(tree) > 1) tree <- paste(tree, collapse = "")

    single_quotes <- function(x, z) {
        x <- charToRaw(x)
        z <- which(x == as.raw(39))
        if (length(z) %% 2) stop("wrong number of single quotes around labels")
        l <- length(z) / 2
        opening <- z[c(TRUE, FALSE)]
        closing <- z[c(FALSE, TRUE)]
        from <- c(1, closing + 1L)
        to <- c(opening - 1L, length(x))
        i <- mapply(":", from = from, to = to, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        keep <- lapply(i, function(i) x[i])
        tmp_label <- paste0("IMPROBABLEPREFIX", 1:l, "IMPROBABLESUFFIX")
        tmpLabsRaw <- lapply(tmp_label, charToRaw)
        n <- 2 * l + 1L
        res <- vector("list", n)
        res[seq(1, n, 2)] <- keep
        res[seq(2, n - 1, 2)] <- tmpLabsRaw
        tree <<- rawToChar(unlist(res))
        i <- mapply(":", from = opening, to = closing, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        orig_label <- lapply(i, function(i) x[i])
        sapply(orig_label, rawToChar)
    }

    ## replace labels with single quotes (if needed)
    SINGLE.QUOTES.FOUND <- grepl("'", tree)
    if (SINGLE.QUOTES.FOUND) tmp_label <- single_quotes(tree)

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
    STRING <- gsub("^_+|_+$", "", STRING)

    tree <- gsub("[ \t]", "", tree) # spaces and TABs within quoted labels are not deleted

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

    if (SINGLE.QUOTES.FOUND) {
        FOO <- function(x) {
            i <- gsub("^IMPROBABLEPREFIX|IMPROBABLESUFFIX$", "", x)
            tmp_label[as.integer(i)]
        }
        for (i in 1:Ntree) {
            lab <- obj[[i]]$tip.label
            k <- grep("IMPROBABLEPREFIX", lab)
            if (length(k)) {
                lab[k] <- FOO(lab[k])
                obj[[i]]$tip.label <- lab
            }
            lab <- obj[[i]]$node.label
            k <- grep("IMPROBABLEPREFIX", lab)
            if (length(k)) {
                lab[k] <- FOO(lab[k])
                obj[[i]]$node.label <- lab
            }
        }
    }
    if (Ntree == 1 && !keep.multi) obj <- obj[[1]] else {
        if (!is.null(tree.names)) names(obj) <- tree.names
        class(obj) <- "multiPhylo"
    }
    obj
}
