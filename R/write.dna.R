## write.dna.R (2020-11-08)

##   Write DNA Sequences in a File

## Copyright 2003-2020 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

write.dna <- function(x, file, format = "interleaved", append = FALSE,
                      nbcol = 6, colsep = " ", colw = 10, indent = NULL,
                      blocksep = 1)
{
    format <- match.arg(format, c("interleaved", "sequential", "fasta"))
    phylip <- if (format %in% c("interleaved", "sequential")) TRUE else FALSE
    if (inherits(x, "DNAbin")) x <- as.character(x)
    aligned <- TRUE
    if (is.matrix(x)) {
        N <- dim(x)
        S <- N[2]
        N <- N[1]
        xx <- vector("list", N)
        for (i in 1:N) xx[[i]] <- x[i, ]
        names(xx) <- rownames(x)
        x <- xx
        rm(xx)
    } else {
        N <- length(x)
        S <- unique(lengths(x, use.names = FALSE))
        if (length(S) > 1) aligned <- FALSE
    }
    if (is.null(names(x))) names(x) <- as.character(1:N)
    if (is.null(indent))
      indent <- if (phylip) 10 else  0
    if (is.numeric(indent))
        indent <- paste(rep(" ", indent), collapse = "")
    if (format == "interleaved") {
        blocksep <- paste(rep("\n", blocksep), collapse = "")
        if (nbcol < 0) format <- "sequential"
    }
    zz <- if (append) file(file, "a") else file(file, "w")
    on.exit(close(zz))
    if (phylip) {
        if (!aligned)
            stop("sequences must have the same length for
 interleaved or sequential format.")
        cat(N, " ", S, "\n", sep = "", file = zz)
        if (nbcol < 0) {
            nb.block <- 1
            nbcol <- totalcol <- ceiling(S/colw)
        } else {
            nb.block <- ceiling(S/(colw * nbcol))
            totalcol <- ceiling(S/colw)
        }
        ## Prepare the sequences in a matrix whose elements are
        ## strings with `colw' characters.
        SEQ <- matrix("", N, totalcol)
        for (i in 1:N) {
            X <- paste(x[[i]], collapse = "")
            for (j in 1:totalcol)
                SEQ[i, j] <- substr(X, 1 + (j - 1)*colw, colw + (j - 1)*colw)
        }
        ## Prepare the names so that they all have the same nb of chars
        max.nc <- max(nchar(names(x)))
        ## always put a space between the sequences and the taxa names
        fmt <- paste("%-", max.nc + 1, "s", sep = "")
        names(x) <- sprintf(fmt, names(x))
    }
    switch(format, "interleaved" = {
        ## Write the first block with the taxon names
        colsel <- if (nb.block == 1) 1:totalcol else 1:nbcol
        for (i in 1:N) {
            cat(names(x)[i], file = zz)
            cat(SEQ[i, colsel], sep = colsep, file = zz)
            cat("\n", file = zz)
        }
        ## Write eventually the other blocks
        if (nb.block > 1) {
            for (k in 2:nb.block) {
                cat(blocksep, file = zz)
                endcolsel <- if (k == nb.block) totalcol else nbcol + (k - 1)*nbcol
                for (i in 1:N) {
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1)*nbcol):endcolsel], sep = colsep, file = zz)
                    cat("\n", file = zz)
                }
            }
        }

    }, "sequential" = {
        if (nb.block == 1) {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, ], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        } else {
            for (i in 1:N) {
                cat(names(x)[i], file = zz)
                cat(SEQ[i, 1:nbcol], sep = colsep, file = zz)
                cat("\n", file = zz)
                for (k in 2:nb.block) {
                    endcolsel <- if (k == nb.block) totalcol else nbcol + (k - 1)*nbcol
                    cat(indent, file = zz)
                    cat(SEQ[i, (1 + (k - 1)*nbcol):endcolsel], sep = colsep, file = zz)
                    cat("\n", file = zz)
                }
            }
        }
    }, "fasta" = {
        for (i in 1:N) {
            cat(">", names(x)[i], file = zz, sep = "")
            cat("\n", file = zz)
            X <- paste(x[[i]], collapse = "")
            S <- length(x[[i]])
            totalcol <- ceiling(S/colw)
            if (nbcol < 0) nbcol <- totalcol
            nb.lines <- ceiling(totalcol/nbcol)
            SEQ <- character(totalcol)
            for (j in 1:totalcol)
                SEQ[j] <- substr(X, 1 + (j - 1) * colw, colw + (j - 1) * colw)
            for (k in 1:nb.lines) {
                endsel <-
                    if (k == nb.lines) length(SEQ) else nbcol + (k - 1)*nbcol
                cat(indent, file = zz)
                cat(SEQ[(1 + (k - 1)*nbcol):endsel], sep = colsep, file = zz)
                cat("\n", file = zz)
            }
        }
    })
}

write.FASTA <- function(x, file, header = NULL, append = FALSE)
{
    dna <- inherits(x, "DNAbin")
    if (!dna && !inherits(x, "AAbin"))
        stop("data are apparently neither DNA nor AA sequences")
    if (!is.null(header)) {
        header <- as.character(header)
        if (!length(header) || !sum(nchar(header)) || is.na(header)) {
            warning("header cannot be coerced as character; was ignored")
            header <- NULL
        }
    }
    labs <- labels(x)
    if (is.matrix(x)) {
        s <- ncol(x) # always integer
        n <- nrow(x) #
    } else {
        s <- -1L
        n <- length(x)
    }
    if (is.null(labs)) labs <- as.character(1:n)
    labs <- lapply(labs, charToRaw)
    if (!is.null(header)) {
        cat(header, sep = "\n", file = file, append = append)
    } else {
        fileExists <- file.exists(file)
        if (append) {
            if (!fileExists) stop(paste("cannot access FASTA file", file))
        } else {
            if (fileExists) file.remove(file)
            o <- file.create(file)
            if (!o) stop(paste("cannot create FASTA file", file))
        }
    }
    ## 'file' should always exist now
    file <- normalizePath(file)
    if (dna) .Call(writeDNAbinToFASTA, x, file, n, s, labs)
    else .Call(writeAAbinToFASTA, x, file, n, s, labs)
    invisible(NULL)
}
