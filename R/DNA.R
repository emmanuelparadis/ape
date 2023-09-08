## DNA.R (2023-02-13)

##   Manipulations and Comparisons of DNA and AA Sequences

## Copyright 2002-2023 Emmanuel Paradis, 2015 Klaus Schliep, 2017 Franz Krah

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

DNAbin2indel <- function(x)
{
    if (is.list(x)) x <- as.matrix(x)
    d <- dim(x)
    s <- as.integer(d[2])
    n <- as.integer(d[1])
    if (s * n > 2^31 - 1)
        stop("DNAbin2indel() cannot handle more than 2^31 - 1 bases")
    res <- .C(DNAbin2indelblock, x, n, s, integer(n*s), NAOK = TRUE)[[4]]
    dim(res) <- d
    rownames(res) <- rownames(x)
    res
}

labels.DNAbin <- function(object, ...)
{
    if (is.list(object)) return(names(object))
    if (is.matrix(object)) return(rownames(object))
    NULL
}

del.gaps <- function(x)
{
    deleteGaps <- function(x) {
        i <- which(x == 4)
        if (length(i)) x[-i] else x
    }

    if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
    if (is.matrix(x)) {
        n <- dim(x)[1]
        y <- vector("list", n)
        for (i in 1:n) y[[i]] <- x[i, ]
        names(y) <- rownames(x)
        x <- y
        rm(y)
    }
    if (!is.list(x)) return(deleteGaps(x))
    x <- lapply(x, deleteGaps)
    class(x) <- "DNAbin"
    x
}

del.rowgapsonly <- function(x, threshold = 1, freq.only = FALSE)
{
    if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
    if (!is.matrix(x)) stop("DNA sequences not in a matrix")
    foo <- function(x) sum(x == 4)
    g <- apply(x, 1, foo)
    if (freq.only) return(g)
    i <- which(g / ncol(x) >= threshold)
    if (length(i)) x <- x[-i, ]
    x
}

del.colgapsonly <- function(x, threshold = 1, freq.only = FALSE)
{
    if (!inherits(x, "DNAbin")) x <- as.DNAbin(x)
    if (!is.matrix(x)) stop("DNA sequences not in a matrix")
    foo <- function(x) sum(x == 4)
    g <- apply(x, 2, foo)
    if (freq.only) return(g)
    i <- which(g / nrow(x) >= threshold)
    if (length(i)) x <- x[, -i]
    x
}

as.alignment <- function(x)
{
    if (is.list(x)) n <- length(x)
    if (is.matrix(x)) n <- dim(x)[1]
    seq <- character(n)
    if (is.list(x)) {
        nam <- names(x)
        for (i in 1:n)
          seq[i] <- paste(x[[i]], collapse = "")
    }
    if (is.matrix(x)) {
        nam <- dimnames(x)[[1]]
        for (i in 1:n)
          seq[i] <- paste(x[i, ], collapse = "")
    }
    obj <- list(nb = n, seq = seq, nam = nam, com = NA)
    class(obj) <- "alignment"
    obj
}

"[.DNAbin" <- function(x, i, j, drop = FALSE)
{
    ans <- NextMethod("[", drop = drop)
    class(ans) <- "DNAbin"
    ans
}

as.matrix.DNAbin <- function(x, ...)
{
    if (is.matrix(x)) return(x)
    if (!is.list(x)) { # vector
        dim(x) <- c(1, length(x))
        return(x)
    }
    s <- unique(lengths(x, use.names = FALSE))
    if (length(s) != 1)
        stop("DNA sequences in list not of the same length.")
    n <- length(x)
    y <- matrix(raw(), n, s)
    for (i in seq_len(n)) y[i, ] <- x[[i]]
    rownames(y) <- names(x)
    class(y) <- "DNAbin"
    y
}

as.list.DNAbin <- function(x, ...)
{
    if (is.list(x)) return(x)
    if (is.null(dim(x))) obj <- list(x) # cause is.vector() doesn't work
    else { # matrix
        class(x) <- NULL
        n <- nrow(x)
        obj <- vector("list", n)
        for (i in seq_len(n)) obj[[i]] <- x[i, , drop = TRUE]
        names(obj) <- rownames(x)
    }
    class(obj) <- "DNAbin"
    obj
}

rbind.DNAbin <- function(...)
{
    obj <- list(...)
    n <- length(obj)
    if (n == 1) return(obj[[1]])
    for (i in 1:n)
        if (!is.matrix(obj[[1]]))
            stop("the 'rbind' method for \"DNAbin\" accepts only matrices")
    NC <- unlist(lapply(obj, ncol))
    if (length(unique(NC)) > 1)
        stop("matrices do not have the same number of columns.")
    for (i in 1:n) class(obj[[i]]) <- NULL # safe but maybe not really needed
    structure(do.call(rbind, obj), class = "DNAbin")
}

cbind.DNAbin <-
    function(..., check.names = TRUE, fill.with.gaps = FALSE,
             quiet = FALSE)
{
    obj <- list(...)
    n <- length(obj)
    if (n == 1) return(obj[[1]])
    for (i in 1:n)
        if (!is.matrix(obj[[1]]))
            stop("the 'cbind' method for \"DNAbin\" accepts only matrices")
    NR <- unlist(lapply(obj, nrow))
    for (i in 1:n) class(obj[[i]]) <- NULL
    if (check.names) {
        NMS <- lapply(obj, rownames)
        for (i in 1:n)
            if (anyDuplicated(NMS[[i]]))
                stop("Duplicated rownames in matrix ", i, ": see ?cbind.DNAbin")
        nms <- unlist(NMS)
        if (fill.with.gaps) {
            NC <- unlist(lapply(obj, ncol))
            nms <- unique(nms)
            ans <- matrix(as.raw(4), length(nms), sum(NC))
            rownames(ans) <- nms
            from <- 1
            for (i in 1:n) {
                to <- from + NC[i] - 1
                k <- match(NMS[[i]], nms)
                ans[k, from:to] <- obj[[i]]
                from <- to + 1
            }
        } else {
            tab <- table(nms)
            ubi <- tab == n
            nms <- names(tab)[which(ubi)]
            ans <- obj[[1]][nms, , drop = FALSE]
            for (i in 2:n)
                ans <- cbind(ans, obj[[i]][nms, , drop = FALSE])
            if (!quiet && !all(ubi))
                warning("some rows were dropped.")
        }
    } else {
        if (length(unique(NR)) > 1)
            stop("matrices do not have the same number of rows.")
        ans <- matrix(unlist(obj), NR)
        rownames(ans) <- rownames(obj[[1]])
    }
    class(ans) <- "DNAbin"
    ans
}

c.DNAbin <- function(..., recursive = FALSE)
{
    if (!all(unlist(lapply(list(...), is.list))))
        stop("the 'c' method for \"DNAbin\" accepts only lists")
    structure(NextMethod("c"), class = "DNAbin")
}

print.DNAbin <- function(x, printlen = 6, digits = 3, ...)
{
    if (is.list(x)) {
        n <- length(x)
        nms <- names(x)
        if (n == 1) {
            cat("1 DNA sequence in binary format stored in a list.\n\n")
            nTot <- length(x[[1]])
            cat("Sequence length:", nTot, "\n")
        } else {
            cat(n, "DNA sequences in binary format stored in a list.\n\n")
            tmp <- lengths(x, use.names = FALSE)
            nTot <- sum(as.numeric(tmp))
            mini <- min(tmp)
            maxi <- max(tmp)
            if (mini == maxi)
                cat("All sequences of same length:", maxi, "\n")
            else {
                cat("Mean sequence length:", round(mean(tmp), 3), "\n")
                cat("   Shortest sequence:", mini, "\n")
                cat("    Longest sequence:", maxi, "\n")
            }
        }
    } else {
        nTot <- length(x)
        if (is.matrix(x)) {
            nd <- dim(x)
            n <- nd[1]
            nms <- rownames(x)
            if (n == 1) {
                cat("1 DNA sequence in binary format stored in a matrix.\n\n")
                cat("Sequence length:", nd[2], "\n")
            } else {
                cat(n, "DNA sequences in binary format stored in a matrix.\n\n")
                cat("All sequences of same length:", nd[2], "\n")
            }
        } else {
            cat("1 DNA sequence in binary format stored in a vector.\n\n")
            cat("Sequence length:", nTot, "\n\n")
        }
    }
    if (exists("nms")) {
        HEAD <- if (n == 1) "\nLabel:" else "\nLabels:"
        TAIL <- ""
        if (printlen < n) {
            nms <- nms[1:printlen]
            TAIL <- "...\n"
        }
        if (any(longs <- nchar(nms) > 60))
            nms[longs] <- paste0(substr(nms[longs], 1, 60), "...")
        cat(HEAD, nms, TAIL, sep = "\n")
    }
    if (nTot <= 1e7) {
        cat("Base composition:\n")
        print(round(base.freq(x), digits))
    } else {
        cat("More than 10 million bases: not printing base composition.\n")
    }
    if (nTot > 1) {
        k <- floor(log(nTot, 1000))
        units <- c("bases", "kb", "Mb", "Gb", "Tb", "Pb", "Eb")
        cat("(Total: ", round(nTot/1000^k, 2), " ", units[k + 1], ")\n", sep = "")
    }
}

as.DNAbin <- function(x, ...) UseMethod("as.DNAbin")

._cs_ <- c("a", "g", "c", "t", "r", "m", "w", "s", "k",
           "y", "v", "h", "d", "b", "n", "-", "?")

._bs_ <- c(136, 72, 40, 24, 192, 160, 144, 96, 80,
           48, 224, 176, 208, 112, 240, 4, 2)

## by Klaus:
as.DNAbin.character <- function(x, ...)
{
    ans <- as.raw(._bs_)[match(tolower(x), ._cs_)]
    if (is.matrix(x)) {
        dim(ans) <- dim(x)
        dimnames(ans) <- dimnames(x)
    }
    class(ans) <- "DNAbin"
    ans
}

as.DNAbin.alignment <- function(x, ...)
{
    n <- x$nb
    x$seq <- tolower(x$seq)
    ans <- matrix("", n, nchar(x$seq[1]))
    for (i in 1:n)
        ans[i, ] <- strsplit(x$seq[i], "")[[1]]
    rownames(ans) <- gsub(" +$", "", gsub("^ +", "", x$nam))
    as.DNAbin.character(ans)
}

as.DNAbin.list <- function(x, ...)
{
    obj <- lapply(x, as.DNAbin)
    class(obj) <- "DNAbin"
    obj
}

as.character.DNAbin <- function(x, ...)
{
    f <- function(xx) {
        ans <- ._cs_[match(as.numeric(xx), ._bs_)]
        if (is.matrix(xx)) {
            dim(ans) <- dim(xx)
            dimnames(ans) <- dimnames(xx)
        }
        ans
    }
    if (is.list(x)) lapply(x, f) else f(x)
}

base.freq <- function(x, freq = FALSE, all = FALSE)
{
    if (!inherits(x, "DNAbin"))
        stop('base.freq requires an object of class "DNAbin"')

    f <- function(x) .Call(BaseProportion, x)

    if (is.list(x)) {
        BF <- rowSums(sapply(x, f))
        n <- sum(as.double(lengths(x, use.names = FALSE)))
    } else {
        n <- length(x)
        BF <- f(x)
    }

    names(BF) <- c("a", "c", "g", "t", "r", "m", "w", "s",
                   "k", "y", "v", "h", "d", "b", "n", "-", "?")
    if (all) {
        if (!freq) BF <- BF / n
    } else {
        BF <- BF[1:4]
        if (!freq) BF <- BF / sum(BF)
    }
    BF
}

Ftab <- function(x, y = NULL)
{
    if (is.null(y)) {
        if (is.list(x)) {
            y <- x[[2]]
            x <- x[[1]]
            if (length(x) != length(y))
                stop("'x' and 'y' not of the same length")
        } else { # 'x' is a matrix
            y <- x[2, , drop = TRUE]
            x <- x[1, , drop = TRUE]
        }
    } else {
        x <- as.vector(x)
        y <- as.vector(y)
        if (length(x) != length(y))
            stop("'x' and 'y' not of the same length")
    }
    out <- matrix(0, 4, 4)
    k <- c(136, 40, 72, 24)
    for (i in 1:4) {
        a <- x == k[i]
        for (j in 1:4) {
            b <- y == k[j]
            out[i, j] <- sum(a & b)
        }
    }
    dimnames(out)[1:2] <- list(c("a", "c", "g", "t"))
    out
}

GC.content <- function(x) sum(base.freq(x)[2:3])

seg.sites <- function(x, strict = FALSE, trailingGapsAsN = TRUE)
{
    if (is.list(x)) x <- as.matrix(x)
    ## is.vector() returns FALSE because of the class,
    ## so we use a different test
    dx <- dim(x)
    if (is.null(dx)) return(integer())
    if (dx[1] == 1) return(integer())
    if (trailingGapsAsN) x <- latag2n(x)
    ans <- .Call(SegSites, x, strict)
    which(as.logical(ans))
}

dist.dna <- function(x, model = "K80", variance = FALSE, gamma = FALSE,
                     pairwise.deletion = FALSE, base.freq = NULL,
                     as.matrix = FALSE)
{
    MODELS <- c("RAW", "JC69", "K80", "F81", "K81", "F84", "T92", "TN93",
                "GG95", "LOGDET", "BH87", "PARALIN", "N", "TS", "TV",
                "INDEL", "INDELBLOCK")
    imod <- pmatch(toupper(model), MODELS)
    if (is.na(imod))
        stop(paste("'model' must be one of:",
                   paste("\"", MODELS, "\"", sep = "", collapse = " ")))
    if (imod == 11 && variance) {
        warning("computing variance not available for model BH87")
        variance <- FALSE
    }
    if (gamma && imod %in% c(1, 5:7, 9:17)) {
        warning(paste("gamma-correction not available for model", model))
        gamma <- FALSE
    }
    if (is.list(x)) x <- as.matrix(x)
    nms <- dimnames(x)[[1]]
    n <- dim(x)[1] # in case nms is NULL

    if (imod %in% c(4, 6:8)) {
        BF <- if (is.null(base.freq)) base.freq(x) else base.freq
    } else BF <- 0

    if (imod %in% 16:17) pairwise.deletion <- TRUE

    if (!pairwise.deletion) {
        keep <- .Call(GlobalDeletionDNA, x)
        x <- x[, as.logical(keep)]
    }
    if (!gamma) {
        alpha <- 0
    } else {
        alpha <- gamma
        gamma <- 1L
    }
    d <- .Call(dist_dna, x, imod, BF, as.integer(pairwise.deletion),
               as.integer(variance), as.integer(gamma), alpha)
    if (variance) {
        var <- d[[2]]
        d <- d[[1]]
    }
    if (imod == 11) {
        dim(d) <- c(n, n)
        dimnames(d) <- list(nms, nms)
    } else {
        attr(d, "Size") <- n
        attr(d, "Labels") <- nms
        attr(d, "Diag") <- attr(d, "Upper") <- FALSE
        attr(d, "call") <- match.call()
        attr(d, "method") <- model
        class(d) <- "dist"
        if (as.matrix) d <- as.matrix(d)
    }
    if (variance) attr(d, "variance") <- var
    d
}

alview <- function(x, file = "", uppercase = TRUE, showpos = TRUE)
{
    if (is.list(x)) x <- as.matrix(x)
    taxa <- formatC(labels(x), width = -1)
    x <- as.character(x)
    s <- ncol(x)
    if (nrow(x) > 1) {
        for (j in seq_len(s)) {
            q <- which(x[-1L, j] == x[1L, j]) + 1L
            x[q, j] <- "."
        }
    }
    x <- apply(x, 1L, paste, collapse = "")
    if (uppercase) x <- toupper(x)
    res <- paste(taxa, x)
    if ((is.logical(showpos) && showpos) || is.numeric(showpos)) {
        if (is.logical(showpos)) {
            pos <- 1:s
            digits <- floor(log10(s)) + 1
        } else {
            pos <- showpos
            digits <- floor(log10(max(pos))) + 1
        }
        hdr <- sprintf(paste0("%0", digits, "d"), pos)
        hdr <- unlist(strsplit(hdr, ""))
        dim(hdr) <- c(digits, length(pos))
        hdr <- apply(hdr, 1, paste, collapse = "")
        hdr <- formatC(hdr, width = nchar(res[1]))
        cat(hdr, file = file, sep = "\n")
    }
    cat(res, file = file, sep = "\n", append = TRUE)
}

where <- function(x, pattern)
{
    pat <- strsplit(pattern, NULL)[[1]]
    if (inherits(x, "DNAbin")) {
        pat <- as.DNAbin(pat)
    } else {
        if (inherits(x, "AAbin")) {
            pat <- as.AAbin(toupper(pat))
        } else {
            stop("'x' should inherit class \"DNAbin\" or \"AAbin\"")
        }
    }
    p <- length(pat)

    f <- function(x, pat, p) {
        if (length(x) < p) {
            warning("sequence shorter than the pattern: returning NULL")
            return(NULL)
        }
        .Call(C_where, x, pat)
    }

    if (is.list(x)) return(lapply(x, f, pat = pat, p = p))
    if (is.matrix(x)) {
        n <- nrow(x)
        res <- vector("list", n)
        for (i in seq_len(n))
            res[[i]] <- f(x[i, , drop = TRUE], pat, p)
        names(res) <- rownames(x)
        return(res)
    }
    f(x, pat, p) # if x is a vector
}

## conversions from BioConductor:

## DNA:
.DNAString2DNAbin <- function(from)
    .Call("charVectorToDNAbinVector", as.character(from))

as.DNAbin.DNAString <- function(x, ...)
{
    res <- list(.DNAString2DNAbin(x))
    class(res) <- "DNAbin"
    res
}

as.DNAbin.DNAStringSet <- function(x, ...)
{
    res <- lapply(x, .DNAString2DNAbin)
    class(res) <- "DNAbin"
    res
}

as.DNAbin.DNAMultipleAlignment <- function(x, ...)
    as.matrix(as.DNAbin.DNAStringSet(as(x, "DNAStringSet")))

as.DNAbin.PairwiseAlignmentsSingleSubject <- function(x, ...)
    as.DNAbin.DNAMultipleAlignment(x)

## AA:
.AAString2AAbin <- function(from)
    charToRaw(as.character(from))

as.AAbin.AAString <- function(x, ...)
{
    res <- list(.AAString2AAbin(x))
    class(res) <- "AAbin"
    res
}

as.AAbin.AAStringSet <- function(x, ...)
{
    res <- lapply(x, .AAString2AAbin)
    class(res) <- "AAbin"
    res
}

as.AAbin.AAMultipleAlignment <- function(x, ...)
    as.matrix(as.AAbin.AAStringSet(as(x, "AAStringSet")))

complement <- function(x)
{
    f <- function(x) {
        ## reorder the vector of raws to match the complement:
        comp <- as.raw(._bs_[c(4:1, 10:9, 7:8, 6:5, 14:11, 15:17)])
        ans <- comp[match(as.integer(x), ._bs_)]
        rev(ans) # reverse before returning
    }
    if (is.matrix(x)) {
        for (i in 1:nrow(x)) x[i, ] <- f(x[i, ])
        return(x)
    } else if (is.list(x)) {
        x <- lapply(x, f)
    } else x <- f(x)
    class(x) <- "DNAbin"
    x
}

trans <- function(x, code = 1, codonstart = 1)
{
    f <- function(x, s, code)
        .C(trans_DNA2AA, x, as.integer(s), raw(s/3), as.integer(code),
           NAOK = TRUE)[[3]]
    if (code > 6)
        stop("only the genetic codes 1--6 are available for now")

    if (codonstart > 1) {
        del <- -(1:(codonstart - 1))
        if (is.list(x)) {
            for (i in seq_along(x)) x[[i]] <- x[[i]][del]
        } else {
            x <- if (is.matrix(x)) x[, del] else x[del]
        }
    }

    if (is.list(x)) {
        res <- lapply(x, trans, code = code)
    } else {
        s <- if (is.matrix(x)) ncol(x) else length(x)
        rest <- s %% 3
        if (rest != 0) {
            s <- s - rest
            x <- if (is.matrix(x)) x[, 1:s] else x[1:s]
            msg <- paste("sequence length not a multiple of 3:", rest, "nucleotide")
            if (rest == 2) msg <- paste0(msg, "s")
            warning(paste(msg, "dropped"))
        }
        if (is.matrix(x)) {
            res <- t(apply(x, 1, f, s = s, code = code))
            if (s == 3) {
                res <- t(res)
                rownames(res) <- rownames(x)
            }
        } else {
            res <- f(x, s, code)
        }
    }
    class(res) <- "AAbin"
    res
}

print.AAbin <- function(x, ...)
{
    if (is.list(x)) {
        n <- length(x)
        cat(n, "amino acid sequence")
        if (n > 1) cat("s")
        cat(" in a list\n\n")
        tmp <- lengths(x, use.names = FALSE)
        maxi <- max(tmp)
        mini <- min(tmp)
        if (mini == maxi)
            cat("All sequences of the same length:", maxi, "\n")
        else {
            cat("Mean sequence length:", round(mean(tmp), 3),
                "\n   Shortest sequence:", mini,
                "\n    Longest sequence:", maxi, "\n")
        }
    } else if (is.matrix(x)) {
        n <- nrow(x)
        cat(n, "amino acid sequence")
        if (n > 1) cat("s")
        cat(" in a matrix\n")
        if (n == 1) cat("Sequence length: ")
        else cat("All sequences of the same length: ")
        cat(ncol(x), "\n")
    } else {
        cat("1 amino acid sequence in a vector:\n\n",
            rawToChar(x))
    }
    cat("\n")
}

"[.AAbin" <- function (x, i, j, drop = FALSE)
{
    ans <- NextMethod("[", drop = drop)
    class(ans) <- "AAbin"
    ans
}

as.character.AAbin <- function(x, ...)
{
    f <- function(xx) {
        ans <- strsplit(rawToChar(xx), "")[[1]]
        if (is.matrix(xx)) {
            dim(ans) <- dim(xx)
            dimnames(ans) <- dimnames(xx)
        }
        ans
    }
    if (is.list(x)) lapply(x, f) else f(x)
}

as.AAbin <- function(x, ...) UseMethod("as.AAbin")

as.AAbin.character <- function(x, ...)
{
    f <- function(x) charToRaw(paste(x, collapse = ""))
    res <- if (is.vector(x)) f(x) else t(apply(x, 1, f))
    class(res) <- "AAbin"
    res
}

labels.AAbin <- function(object, ...) labels.DNAbin(object, ...)

## TO BE MOVED TO phangorn LATER
if (getRversion() >= "2.15.1") utils::globalVariables("phyDat")
as.phyDat.AAbin <- function(x, ...) phyDat(as.character(x), type = "AA")
## \alias{as.phyDat.AAbin}
## \method{as.phyDat}{AAbin}(x, \dots)

dist.aa <- function(x, pairwise.deletion = FALSE, scaled = FALSE)
{
    n <- nrow(x)
    d <- numeric(n*(n - 1)/2)
    X <- charToRaw("X")
    k <- 0L
    if (!pairwise.deletion) {
        del <- apply(x, 2, function(y) any(y == X))
        if (any(del)) x <- x[, !del]
        for (i in 1:(n - 1)) {
            for (j in (i + 1):n) {
                k <- k + 1L
                d[k] <- sum(x[i, ] != x[j, ])
            }
        }
        if (scaled) d <- d/ncol(x)
    } else {
        for (i in 1:(n - 1)) {
            a <- x[i, ]
            for (j in (i + 1):n) {
                b <- x[j, ]
                del <- a == X | b == X
                p <- length(b <- b[!del])
                tmp <- sum(a[!del] != b)
                k <- k + 1L
                d[k] <- if (scaled) tmp/p else tmp
            }
        }
    }
    attr(d, "Size") <- n
    attr(d, "Labels") <- rownames(x)
    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    d
}

AAsubst <- function(x)
{
    X <- charToRaw("X")
    f <- function(y) length(unique.default(y[y != X]))
    which(apply(x, 2, f) > 1)
}

.AA_3letter <- c("Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile",  "Lys", "Leu",
                 "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr",  "Val", "Trp", "Tyr",
                 "Xaa", "Stp")

.AA_1letter <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                 "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "X", "*")

.AA_raw <- sapply(.AA_1letter, charToRaw)

.AA_3cat <- list(Hydrophobic = .AA_raw[c("V", "I", "L", "F", "W", "Y", "M")],
                 Small = .AA_raw[c("P", "G", "A", "C")],
                 Hydrophilic = .AA_raw[c("S", "T", "H", "N", "Q", "D", "E", "K", "R")])

#.AA_3cat <- list(Hydrophobic = .AA_raw[c("V", "I", "L", "F", "W", "Y", "M")],
#                 Small = .AA_raw[c("P", "G", "A", "C")],
#                 Hydrophilic = .AA_raw[c("S", "T", "H", "N", "Q", "D", "E", "K", "R")])



Ape_NT <- list(properties = list(
    a="a", g="g", c="c", t="t", n="n", "-"="-"),
    color=c("red", "yellow", "green", "blue", "grey", "black"))


RY_NT <- list(properties = list(
    Purine = c("a", "g", "r"),
    Pyrimidine = c("c", "t", "y"),
    "n" = "n",
    "-" = "-"),
    color=c("#FF00FF", "#00FFFF", "grey", "black"))



Ape_AA <- list(properties = list(
    Hydrophobic = c("V", "I", "L", "F", "W", "Y", "M"),
    Small = c("P", "G", "A", "C"),
    Hydrophilic = c("S", "T", "H", "N", "Q", "D", "E", "K", "R")),
    color=c("red", "yellow", "blue"))

# Properties + Conservation (Clustal X)
Clustal <- list(properties = list(
    Hydrophobic = c("A", "I", "L", "M", "F", "W", "V"),
    Positive = c("K", "R"),
    Negative = c("E", "D"),
    Polar = c("N", "Q", "S", "T"),
    Glycines = "G",
    Prolines = "P",
    Aromatic = c("H", "Y"),
    Cysteine = "C"),
    color= c("#80a0f0", "#f01505", "#c048c0", "#15c015", "#f09048", "#c0c000",
             "#15a4a4", "#f08080")
)


# Polarity geneious
Polarity <- list(properties = list(
    "Non polar" = c("G", "A", "V", "L", "I", "F", "W", "M", "P"),
    "Polar, uncharged" = c("S", "T", "C", "Y", "N", "Q"),
    "Polar, acidic" = c("D", "E"),
    "Polar, basic" = c("K", "R", "H")),
    color = c("yellow", "green", "red", "blue"))

# Physicochemical Properties
Zappo_AA <- list(properties = list(
    Hydrophobic = c("I", "L", "V", "A", "M"), # `Aliphatic/Hydrophobic`
    Aromatic = c("F", "W", "Y"),
    Positive = c("K", "R", "H"),
    Negative = c("E", "D"),
    Hydrophilic = c("S", "T", "N", "Q"),
    Conformational = c("P", "G"), # `Conformationally special`
    Cysteine = "C"),
    color= c("#ff7979", "#f89f56", "#0070c0", "#c00000", "#08c81a", "#cc00cc",
             "#ffff00")
)

Transmembrane_tendency <- list(properties = list(
    Lys = "K",
    Asp = "D",
    Glu = "E",
    Arg = "R",
    Gln = "Q",
    Asn = "N",
    Pro = "P",
    His = "H",
    Ser = "S",
    Thr = "T",
    Cys = "C",
    Gly = "G",
    Ala = "A",
    Tyr = "Y",
    Met = "M",
    Val = "V",
    Trp = "W",
    Leu = "L",
    Ile = "I",
    Phe = "F"
), color=c("#0000FF", "#0D00F1", "#1A00E4", "#2800D6", "#3500C9", "#4300BB",
           "#5000AE", "#5D00A1", "#6B0093", "#780086", "#860078", "#93006B",
           "#A1005D", "#AE0050", "#BB0043", "#C90035", "#D60028", "#E4001A",
           "#F1000D", "#FF0000"))


image.worker <-
    function(x, what, col, bg = "white", xlab = "", ylab = "",
             show.labels = TRUE, cex.lab = 1, legend = TRUE,
             grid = FALSE, show.bases = FALSE, base.cex = 1,
             base.font = 1, base.col = "black", scheme=NULL, bs=NULL, cs=NULL,...)
{
#    what <-
#        if (missing(what)) c("a", "g", "c", "t", "n", "-") else tolower(what)
    if (missing(col))
            col <- c("red", "yellow", "green", "blue", "grey", "black")
        x <- as.matrix(x) # tests if all sequences have the same length
        n <- (dx <- dim(x))[1] # number of sequences
        s <- dx[2] # number of sites
        y <- integer(N <- length(x))
        ncl <- length(what)
        col <- rep(col, length.out = ncl)
        brks <- 0.5:(ncl + 0.5)
        sm <- 0L
        for (i in ncl:1) {
            k <- bs[match(what[[i]], cs)]
            sel <- which(unclass(x) %in% k)
            if (L <- length(sel)) {
                y[sel] <- i
                sm <- sm + L
            } #else {
#                what <- what[-i]
#                col <- col[-i]
#                brks <- brks[-i]
#            }
        }
        dim(y) <- dx
        ## if there's no 0 in y, must drop 'bg' from the cols passed to image:
        if (sm == N) {
            leg.co <- co <- col
            leg.txt <- toupper(names(what))
        } else {
            co <- c(bg, col)
            leg.txt <- c(toupper(names(what)), "others")
            leg.co <- c(col, bg)
            brks <- c(-0.5, brks)
        }
        yaxt <- if (show.labels) "n" else "s"
        image.default(1:s, 1:n, t(y[n:1, , drop = FALSE]), col = co, xlab = xlab,
                      ylab = ylab, yaxt = yaxt, breaks = brks, ...)
        if (show.labels)
            mtext(rownames(x), side = 2, line = 0.1, at = n:1,
                  cex = cex.lab, adj = 1, las = 1)
        if (legend) {
            psr <- par("usr")
            xx <- psr[2]/2
            yy <- psr[4] * (0.5 + 0.5/par("plt")[4])
            n_col <- length(leg.txt)
            if(n_col > 6) n_col <- ceiling(n_col/2)
            legend(xx, yy, legend = leg.txt, pch = 22, pt.bg = leg.co,
                   pt.cex = 2, bty = "n", xjust = 0.5, yjust = 0.5,
                   horiz = FALSE, ncol=n_col, xpd = TRUE)
        }
        if (grid) {
            if (is.logical(grid)) grid <- 3L
            if (grid %in% 2:3) abline(v = seq(1.5, s - 0.5, 1), lwd = 0.33, xpd = FALSE)
            if (grid %in% c(1, 3)) abline(h = seq(1.5, n - 0.5, 1), lwd = 0.33, xpd = FALSE)
        }
        if (show.bases) {
            x <- toupper(as.character(x))
            xx <- rep(1:s, each = n)
            yy <- rep(n:1, s)
            text(xx, yy, x, cex = base.cex, font = base.font, col = base.col)
        }
    }

image.AAbin <-
    function(x, what, col, bg = "white", xlab = "", ylab = "",
             show.labels = TRUE, cex.lab = 1, legend = TRUE,
             grid = FALSE, show.aa = FALSE, aa.cex = 1,
             aa.font = 1, aa.col = "black", scheme="Ape_AA", ...)
{
    scheme <- match.arg(scheme, c("Ape_AA", "Clustal", "Zappo_AA",
                                  "Polarity", "Transmembrane_tendency"))
    scheme <- get(scheme, environment(image.AAbin))
    if (missing(what)){
        if(!is.null(scheme)) what <- scheme$properties
        else what <- Ape_AA$properties
    }
    if (missing(col)){
        if(!is.null(scheme)) col <- scheme$color
        else col <- c("red", "yellow", "blue")
    }
    image.worker(x, what, col, bg = bg, xlab = xlab, ylab = ylab,
                     show.labels = show.labels, cex.lab = cex.lab,
                     legend = legend,
                     grid = grid, show.bases = show.aa, base.cex = aa.cex,
                     base.font = aa.font, base.col = aa.col, scheme=scheme,
                     bs=.AA_raw, cs=.AA_1letter,...)
}

image.DNAbin <-
    function(x, what, col, bg = "white", xlab = "", ylab = "",
             show.labels = TRUE, cex.lab = 1, legend = TRUE,
             grid = FALSE, show.bases = FALSE, base.cex = 1,
             base.font = 1, base.col = "black", scheme="Ape_NT", ...)
{
    scheme <- match.arg(scheme, c("Ape_NT", "RY_NT"))
    scheme <- get(scheme, environment(image.AAbin))
    if (missing(what)){
        if(!is.null(scheme)) what <- scheme$properties
        else what <- Ape_AA$properties
    }
    if (missing(col)){
        if(!is.null(scheme)) col <- scheme$color
        else col <- c("red", "yellow", "blue")
    }
    image.worker(x, what, col, bg = bg, xlab = xlab, ylab = ylab,
                 show.labels = show.labels, cex.lab = cex.lab, legend = legend,
                 grid = grid, show.bases = show.bases, base.cex = base.cex,
                 base.font = base.font, base.col = base.col, scheme=scheme,
                 bs=as.raw(._bs_), cs=._cs_,...)
}


checkAlignment <- function(x, check.gaps = TRUE, plot = TRUE, what = 1:4)
{
    cat("\nNumber of sequences:", n <- nrow(x),
        "\nNumber of sites:", s <- ncol(x), "\n")
if (check.gaps) {
    cat("\n")
    y <- DNAbin2indel(x)
    gap.length <- sort(unique.default(y))[-1]
    if (!length(gap.length)) cat("No gap in alignment.\n")
    else {
        rest <- gap.length %% 3
        if (any(cond <- rest > 0)) {
            cat("Some gap lengths are not multiple of 3:", gap.length[cond])
        } else cat("All gap lengths are multiple of 3.")

        tab <- tabulate(y, gap.length[length(gap.length)])
        tab <- tab[gap.length]
        cat("\n\nFrequencies of gap lengths:\n")
        names(tab) <- gap.length
        print(tab)

        ## find gaps on the borders:
        col1 <- unique(y[, 1])
        if (!col1[1]) col1 <- col1[-1]
        if (length(col1))
            cat("   => length of gaps on the left border of the alignment:", unique(col1), "\n")
        else cat("   => no gap on the left border of the alignment\n")
        i <- which(y != 0, useNames = FALSE)
        jcol <- i %/% nrow(y) + 1
        yi <- y[i]
        j <- yi == s - jcol + 1
        if (any(j))
            cat("   => length of gaps on the right border of the alignment:", yi[j], "\n")
        else cat("   => no gap on the right border of the alignment\n")

        ## find base segments:
        A <- B <- numeric()
        for (i in seq_len(n)) {
            j <- which(y[i, ] != 0) # j: start of each gap in the i-th sequence
            if (!length(j)) next
            k <- j + y[i, j] # k: start of each base segment in the i-th sequence
            if (j[1] != 1) k <- c(1, k) else j <- j[-1]
            if (k[length(k)] > s) k <- k[-length(k)] else j <- c(j, s + 1)
            A <- c(A, j)
            B <- c(B, k)
        }
        AB <- unique(cbind(A, B))
        not.multiple.of.3 <- (AB[, 1] - AB[, 2]) %% 3 != 0
        left.border <- AB[, 2] == 1
        right.border <- AB[, 1] == s + 1
        Nnot.mult3 <- sum(not.multiple.of.3)
        cat("\nNumber of unique contiguous base segments defined by gaps:", nrow(AB), "\n")
        if (!Nnot.mult3) cat("All segment lengths multiple of 3.\n")
        else {
            Nleft <- sum(not.multiple.of.3 & left.border)
            Nright <- sum(not.multiple.of.3 & right.border)
            cat("Number of segment lengths not multiple of 3:", Nnot.mult3, "\n",
                "   => on the left border of the alignement:", Nleft, "\n",
                "   => on the right border                 :", Nright, "\n")
            if (Nright + Nleft < Nnot.mult3) {
                cat("    => positions of these segments inside the alignment: ")
                sel <- not.multiple.of.3 & !left.border & !right.border
                cat(paste(AB[sel, 2], AB[sel, 1] - 1, sep = ".."), "\n")
            }
        }
    }
} else gap.length <- numeric()
    ss <- seg.sites(x)
    cat("\nNumber of segregating sites (including gaps):", length(ss))
    BF.col <- matrix(NA_real_, length(ss), 4)
    for (i in seq_along(ss))
        BF.col[i, ] <- base.freq(x[, ss[i]])#, freq = TRUE)
    tmp <- apply(BF.col, 1, function(x) sum(x > 0))
    cat("\nNumber of sites with at least one substitution:", sum(tmp > 1))
    cat("\nNumber of sites with 1, 2, 3 or 4 observed bases:\n")
    tab2 <- tabulate(tmp, 4L)
    tab2[1] <- s - sum(tab2)
    names(tab2) <- 1:4
    print(tab2)
    cat("\n")
    H <- numeric(s)
    H[ss] <- apply(BF.col, 1, function(x) {x <- x[x > 0]; -sum(x * log(x))})

    G <- rep(1, s)
    G[ss] <- tmp
    if (plot) {
        if (length(what) == 4) {
            mat <- if (length(gap.length)) 1:4 else c(1, 0, 2, 3)
            layout(matrix(mat, 2, 2))
        } else {
            if (length(what) != 1) {
                what <- what[1]
                warning("argument 'what' has length > 1: the first value is taken")
            }
        }
        if (1 %in% what) image(x)
        if (2 %in% what && length(gap.length)) barplot(tab, xlab = "Gap length")
        if (3 %in% what)
            plot(1:s, H, "h", xlab = "Sequence position", ylab = "Shannon index (H)")
        if (4 %in% what)
            plot(1:s, G, "h", xlab = "Sequence position", ylab = "Number of observed bases")
    }
}

all.equal.DNAbin <- function(target, current, plot = FALSE, ...)
{
    if (identical(target, current)) return(TRUE)
    name.target <- deparse(substitute(target))
    name.current <- deparse(substitute(current))
    st1 <- "convert list as matrix for further comparison."
#    st2 <- ""
    st3 <- "Subset your data for further comparison."
    isali1 <- is.matrix(target)
    isali2 <- is.matrix(current)
    if (isali1 && !isali2)
        return(c("1st object is a matrix, 2nd object is a list:", st1))
    if (!isali1 && isali2)
        return(c("1st object is a list, 2nd object is a matrix:", st1))
    if (!isali1 && !isali2)
        return(c("Both objects are lists:",
                 "convert them as matrices for further comparison."))
#    n1 <- if (isali1) nrow(target) else length(target)
#    n2 <- if (isali2) nrow(current) else length(current)
    if (ncol(target) != ncol(current))
        return("Numbers of columns different: comparison stopped here.")

    foo <- function(n) ifelse(n == 1, "sequence", "sequences")
    doComparison <- function(target, current)
        which(target != current, arr.ind = TRUE, useNames = FALSE)

    n1 <- nrow(target)
    n2 <- nrow(current)
    labs1 <- labels(target)
    labs2 <- labels(current)
    if (identical(labs1, labs2)) {
        res <- "Labels in both objects identical."
        res <- list(messages = res,
                    different.sites = doComparison(target, current))
    } else {
        in12 <- labs1 %in% labs2
        in21 <- labs2 %in% labs1
        if (n1 != n2) {
            res <- c("Number of sequences different:",
                     paste(n1, foo(n1), "in 1st object;", n2,
                           foo(n2), "in 2nd object."),
                     st3)
            plot <- FALSE
        } else { # n1 == n2
            if (any(!in12)) {
                res <- c("X: 1st object (target), Y: 2nd object (current).",
                         paste("labels in X not in Y:",
                               paste(labs1[!in12], collapse = ", ")),
                         paste("labels in X not in Y:",
                               paste(labs2[!in21], collapse = ", ")),
                         st3)
                plot <- FALSE
            } else {
                res <- c("Labels in both objects identical but not in the same order.",
                         "Comparing sequences after reordering rows of the second matrix.")
                current <- current[labs1, ]
                if (identical(target, current)) {
                    res <- c(res, "Sequences are identical.")
                    plot <- FALSE
                } else {
                    res <- list(messages = res,
                                different.sites = doComparison(target, current))
                }
            }
	}
    }
    if (plot) {
        cols <- unique(res$different.sites[, 2])
        diff.cols <- diff(cols)
        j <- which(diff.cols != 1)
        end <- c(cols[j], cols[length(cols)])
        start <- c(cols[1], cols[j + 1])
        v <- cumsum(end - start + 1) + 0.5
        f <- function(lab) {
            axis(2, at = seq_len(n1), labels = FALSE)
            axis(1, at = seq_along(cols), labels = cols)
            mtext(lab, line = 1, adj = 0, font = 2)
        }
        layout(matrix(1:2, 2))
        par(xpd = TRUE)
        image(target[, cols], show.labels = FALSE, axes = FALSE, ...)
        f(name.target)
        xx <- c(0.5, v)
        segments(xx, 0.5, xx, n1, lty = 2, col = "white", lwd = 2)
        segments(xx, 0.5, xx, -1e5, lty = 2, lwd = 2)
        image(current[, cols], show.labels = FALSE, axes = FALSE, ...)
        f(name.current)
        segments(xx, 0.5, xx, n2, lty = 2, col = "white", lwd = 2)
        segments(xx, 1e5, xx, n2, lty = 2, lwd = 2)
        #segments(0.5, -5, length(cols) + 0.5, -5, lwd = 5, col = "grey")
        #rect(0.5, -4, length(cols) + 0.5, -3, col = "grey")
        #segments(0.5, 0.5, 10, -3)
    }
    res
}

## From Franz Krah <f.krah@mailbox.org>:

## estensions of the AAbin class to complement the DNAbin class funcitons

c.AAbin <- function(..., recursive = FALSE)
{
  if (!all(unlist(lapply(list(...), is.list))))
    stop("the 'c' method for \"AAbin\" accepts only lists")
  structure(NextMethod("c"), class = "AAbin")
}

rbind.AAbin <- function(...)
{
    obj <- list(...)
    n <- length(obj)
    if (n == 1) return(obj[[1]])
    for (i in 1:n)
        if (!is.matrix(obj[[1]]))
            stop("the 'rbind' method for \"AAbin\" accepts only matrices")
    NC <- unlist(lapply(obj, ncol))
    if (length(unique(NC)) > 1)
        stop("matrices do not have the same number of columns.")
    for (i in 1:n) class(obj[[i]]) <- NULL # safe but maybe not really needed
    structure(do.call(rbind, obj), class = "AAbin")
}

cbind.AAbin <-
    function(..., check.names = TRUE, fill.with.Xs = FALSE,
             quiet = FALSE)
{
    obj <- list(...)
    n <- length(obj)
    if (n == 1) return(obj[[1]])
    for (i in 1:n)
        if (!is.matrix(obj[[1]]))
            stop("the 'cbind' method for \"AAbin\" accepts only matrices")
    NR <- unlist(lapply(obj, nrow))
    for (i in 1:n) class(obj[[i]]) <- NULL
    if (check.names) {
        NMS <- lapply(obj, rownames)
        for (i in 1:n)
            if (anyDuplicated(NMS[[i]]))
                stop("Duplicated rownames in matrix ", i, ": see ?cbind.AAbin")
        nms <- unlist(NMS)
        if (fill.with.Xs) {
            NC <- unlist(lapply(obj, ncol))
            nms <- unique(nms)
            ans <- matrix(charToRaw("X"), length(nms), sum(NC))
            rownames(ans) <- nms
            from <- 1
            for (i in 1:n) {
                to <- from + NC[i] - 1
                k <- match(NMS[[i]], nms)
                ans[k, from:to] <- obj[[i]]
                from <- to + 1
            }
        } else {
            tab <- table(nms)
            ubi <- tab == n
            nms <- names(tab)[which(ubi)]
            ans <- obj[[1]][nms, , drop = FALSE]
            for (i in 2:n)
                ans <- cbind(ans, obj[[i]][nms, , drop = FALSE])
            if (!quiet && !all(ubi))
                warning("some rows were dropped.")
        }
    } else {
        if (length(unique(NR)) > 1)
            stop("matrices do not have the same number of rows.")
        ans <- matrix(unlist(obj), NR)
        rownames(ans) <- rownames(obj[[1]])
    }
    class(ans) <- "AAbin"
    ans
}

as.AAbin.list <- function(x, ...)
{
  obj <- lapply(x, as.AAbin)
  class(obj) <- "AAbin"
  obj
}

as.list.AAbin <- function(x, ...)
{
  if (is.list(x)) return(x)
  if (is.null(dim(x))) obj <- list(x) # cause is.vector() doesn't work
  else { # matrix
    n <- nrow(x)
    obj <- vector("list", n)
    for (i in seq_len(n)) obj[[i]] <- x[i, , drop = TRUE]
    names(obj) <- rownames(x)
  }
  class(obj) <- "AAbin"
  obj
}

as.matrix.AAbin <- function(x, ...)
{
  if (is.matrix(x)) return(x)
  if (!is.list(x)) { # vector
    dim(x) <- c(1, length(x))
    return(x)
  }
  s <- unique(lengths(x, use.names = FALSE))
  if (length(s) != 1)
    stop("AA sequences in list not of the same length.")
  n <- length(x)
  y <- matrix(raw(), n, s)
  for (i in seq_len(n)) y[i, ] <- x[[i]]
  rownames(y) <- names(x)
  class(y) <- "AAbin"
  y
}

rDNAbin <- function(n, nrow, ncol, base.freq = rep(0.25, 4), prefix = "Ind_")
{
    foo <- function(n, prob) {
        vec <- as.raw(._bs_[1:4])
        vec[sample.int(4L, n, TRUE, prob, FALSE)]
    }
    base.freq <- if (all(base.freq == 0.25)) NULL else base.freq[c(1, 3, 2, 4)]
    if (missing(n)) {
        if (missing(nrow) && missing(ncol))
            stop("nrow and ncol should be given if n is missing")
        res <- foo(nrow * ncol, base.freq)
        dim(res) <- c(nrow, ncol)
        rownames(res) <- paste0(prefix, 1:nrow)
    } else {
        res <- lapply(n, foo, prob = base.freq)
        names(res) <- paste0(prefix, seq_along(n))
    }
    class(res) <- "DNAbin"
    res
}

dnds <- function(x, code = 1, codonstart = 1, quiet = FALSE,
                 details = FALSE, return.categories = FALSE)
{
    if (code > 6) stop("only the genetic codes 1--6 are available for now")
    if (is.list(x)) x <- as.matrix(x)
    n <- nrow(x)
    if (nrow(unique.matrix(x)) != n) stop("sequences are not unique")
### if (any(base.freq(x, TRUE, TRUE)[-(1:4)] > 0)) stop("ambiguous bases are not permitted")
    if (codonstart > 1) {
        del <- -(1:(codonstart - 1))
        x <- x[, del]
    }
    p <- ncol(x)
    rest <- p %% 3
    if (rest) {
        p <- p - rest
        x <- x[, 1:p]
        msg <- sprintf("sequence length not a multiple of 3: %d %s dropped",
                       rest, ngettext(rest, "base", "bases"))
        warning(msg)
    }

    degMat <- .buildDegeneracyMatrix(code)
    Lcat <- matrix(0L, n, p)
    V1 <- V2 <- V3 <- integer(136)
    i <- c(136L, 72L, 40L, 24L)
    V1[i] <- c(1L, 17L, 33L, 49L)
    V2[i] <- c(0L, 4L, 8L, 12L)
    V3[i] <- 0:3

    class(x) <- NULL
    z <- as.integer(x)
    N <- length(x)

    SHIFT <- c(0L, n, 2L * n)
    p <- 1L + SHIFT
    while (p[3] <= N) {
        for (i in 1:n) {
            codon <- z[p]
            ii <- V1[codon[1]] + V2[codon[2]] + V3[codon[3]]
            if (!is.na(ii)) Lcat[p] <- degMat[ii, ]
            p <- p + 1L
        }
        p <- p[3] + SHIFT
    }

    if (return.categories) return(Lcat)
    if (details) quiet <- TRUE

    deg <- c(0, 2, 4) # the 3 levels of degeneracy
    nout <- n*(n - 1)/2
    res <- numeric(nout)
    k <- 1L
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (!quiet) cat("\r", round(100*k/nout), "%")
            z <- x[c(i, j), ]
            Lavg <- (Lcat[i, ] + Lcat[j, ])/2
            Lavg[Lavg == 1] <- 2
            Lavg[Lavg == 3] <- 4
            ii <- lapply(deg, function(x) which(x == Lavg))
            L <- lengths(ii)
            S <- lapply(ii, function(id) dist.dna(z[, id, drop = FALSE], "TS"))
            V <- lapply(ii, function(id) dist.dna(z[, id, drop = FALSE], "TV"))
            S <- unlist(S, use.names = FALSE)
            V <- unlist(V, use.names = FALSE)
            if (details) {
                cat(sprintf("\nComparing sequences %d and %d:\n", i, j))
                tmp <- rbind(S, V)
                dimnames(tmp) <- list(c("Transitions", "Transversions"),
                                      c("Nondegenerate", "Twofold", "Fourfold"))
                print(tmp)
            }
            P <- S/L
            Q <- V/L
            a <- 1/(1 - 2*P - Q)
            b <- 1/(1 - 2*Q)
            c <- (a - b)/2
            A <- log(a)/2 - log(b)/4
            B <- log(b)/2
            dS <- (L[2]*A[2] + L[3]*A[3])/sum(L[2:3]) + B[3]
            dN <- A[1] + (L[1]*B[1] + L[2]*B[2])/sum(L[1:2])
            res[k] <- dN/dS
            k <- k + 1L
        }
    }
    if (!quiet) cat("... done\n")
    attr(res, "Size") <- n
    attr(res, "Labels") <- rownames(x)
    attr(res, "Diag") <- attr(res, "Upper") <- FALSE
    attr(res, "call") <- match.call()
    attr(res, "method") <- "dNdS (Li 1993)"
    class(res) <- "dist"
    res
}

.buildDegeneracyMatrix <- function(code) {
    b <- as.raw(._bs_[1:4])
    CODONS <- cbind(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16))
    AA <- trans(CODONS, code = code)

    degeneracyMatrix <- matrix(0L, 64L, 3L)
    deg <- c(4L, 2L, 2L, 0L)

    ## 1/ find the bases at 3rd positions that are twofold/fourfold degenerate
    s  <- 1:4
    while (s[4L] <= 64) {
        degeneracyMatrix[s, 3L] <- deg[length(unique(AA[s]))]
        s <- s + 4L
    }
    ## 2/ all bases at 2nd positions are nondegenerate: no need to do anything
    ## 3/ are some bases at 1st positions twofold degenerate?
    s <- c(1L, 17L, 33L, 49L)
    while (s[1L] < 17) {
        degeneracyMatrix[s, 1L] <- deg[length(unique(AA[s]))]
        s <- s + 1L
    }
    degeneracyMatrix
}

latag2n <- function(x) {
    if (is.list(x)) x <- as.matrix(x)
    dx <- dim(x)
    clx <- class(x)
    res <- .Call(leading_trailing_gaps_to_N, x)
    ## the order of the next two commands is crucial if 'clx' is
    ## c("matrix", "array") because the dimension must be set before
    ## setting the class (this is the case for pegas::mjn() where the
    ## class "DNAbin" is dropped before calling latag2n())
    dim(res) <- dx
    class(res) <- clx
    res
}

solveAmbiguousBases <- function(x, method = "columnwise", random = TRUE)
{
    if (method == "columnwise") {
        if (is.list(x)) x <- as.matrix(x)
        p <- ncol(x)
        for (j in 1:p) {
            BF <- base.freq(x[, j], TRUE, TRUE)
            ambi <- BF[5:15]
            K <- which(ambi > 0)
            if (length(K)) {
                agct <- BF[c(1, 3, 2, 4)]
                for (b in K) {
                    base <- as.DNAbin(names(ambi[b]))
                    sel <- agct[rev(rawToBits(base))[1:4] == 1]
                    if (!sum(sel)) sel[] <- 1L
                    i <- which(x[, j] == base)
                    tmp <- if (random) sample(names(sel), length(i), TRUE, sel) else names(sel)[which.max(sel)]
                    x[i, j] <- as.DNAbin(tmp)
                }
            }
        }
    }
    x
}

##distK80 <- function(x, pairwise.deletion = FALSE)
##{
##    nms <- dimnames(x)[[1]]
##    n <- length(nms)
##    if (!pairwise.deletion) {
##        keep <- .Call(GlobalDeletionDNA, x)
##        x <- x[, as.logical(keep)]
##        d <- .Call(dist_dna_K80_short, x)
##    } else {
##        d <- .Call(dist_dna_K80_short_pairdel, x)
##    }
##    attr(d, "Size") <- n
##    attr(d, "Labels") <- nms
##    attr(d, "Diag") <- attr(d, "Upper") <- FALSE
##    attr(d, "call") <- match.call()
##    attr(d, "method") <- "K80"
##    class(d) <- "dist"
##    d
##}
