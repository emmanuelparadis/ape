## read.GenBank.R (2020-05-15)

##   Read DNA Sequences and Annotations from GenBank

## Copyright 2002-2020 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

read.GenBank <- function(access.nb, seq.names = access.nb, species.names = TRUE,
                         as.character = FALSE, chunk.size = 400, quiet = TRUE)
{
    chunk.size <- as.integer(chunk.size)
    N <- length(access.nb)
    ## if more than 400 sequences, we break down the requests
    a <- 1L
    b <- if (N > chunk.size) chunk.size else N
    fl <- paste0(tempfile(), ".fas")
    if (!quiet)
        cat("Note: chunk.size =", chunk.size, "(max nb of sequences downloaded together)\n")
    repeat {
        if (!quiet) cat("\rDownloading sequences:", b, "/", N, "...")
        URL <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                      paste(access.nb[a:b], collapse = ","), "&rettype=fasta&retmode=text")
        X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
        cat(X, sep = "\n", file = fl, append = TRUE)
        if (b == N) break
        a <- b + 1L
        b <- b + chunk.size
        if (b > N) b <- N
    }
    if (!quiet) {
        cat(" Done.\nNote: the downloaded sequences are in file:", fl)
        cat("\nReading sequences...")
    }
    res <- read.FASTA(fl)
    if (is.null(res)) return(NULL)
    attr(res, "description") <- names(res)
    if (length(access.nb) != length(res)) {
        names(res) <- gsub("\\..*$", "", names(res))
        failed <- paste(access.nb[! access.nb %in% names(res)], collapse = ", ")
        warning(paste0("cannot get the following sequence(s):\n", failed))
    } else names(res) <- access.nb

    if (as.character) res <- as.character(res)
    if (!quiet) cat("\n")
    if (species.names) {
        a <- 1L
        b <- if (N > chunk.size) chunk.size else N
        sp <- character(0)
        repeat {
            if (!quiet) cat("\rDownloading species names:", b, "/", N)
            URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
                     paste(access.nb[a:b], collapse = ","), "&rettype=gb&retmode=text", sep = "")
            X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE, n = -1)
            sp <- c(sp, gsub(" +ORGANISM +", "", grep("ORGANISM", X, value = TRUE)))
            if (b == N) break
            a <- b + 1L
            b <- b + chunk.size
            if (b > N) b <- N
        }
        if (!quiet) cat(".\n")
        attr(res, "species") <- gsub(" ", "_", sp)
    }
    if (!quiet) cat("Done.\n")
    res
}

.parse.annotations.file <- function(file) {
    get.product.others.gene <- function(a, b) {
        res <- rep(NA_character_, 3L)
        if (a > b) return(res)
        y <- x[a:b]
        li <- length(i <- grep("product\t", y))
        lj <- length(j <- grep("gene\t", y))
        if (li) res[1L] <- gsub("product\t", "", y[i])
        if (lj) res[2L] <- gsub("gene\t", "", y[j])
        if (length(y) > li + lj) res[3L] <- paste(y[-c(i, j)], collapse = "; ")
        res
    }
    convert.with.incomplete <- function(vec) {
        i <- grep("<|>", vec)
        if (length(i)) {
            incomplete <<- c(incomplete, i)
            vec[i] <- gsub("<|>", "", vec[i])
        }
        as.integer(vec)
    }
    x <- scan(file, what = "", sep = "\n", quiet = TRUE, skip = 1)
    n <- length(x)
    i <- grep("^\t\t\t", x)
    i2 <- seq_len(n)[-i]
    Y <- strsplit(x[i2], "\t")

    incomplete <- NULL
    start <- convert.with.incomplete(sapply(Y, "[", 1L))
    end <- convert.with.incomplete(sapply(Y, "[", 2L))
    if (!is.null(incomplete)) {
        incomplete <- sort(unique(incomplete))
        warning(paste("features were incomplete in row(s):",
                      paste(incomplete, collapse = ", ")))
    }

    res <- data.frame(start, end)
    sel <- which(!duplicated.data.frame(res))
    res <- res[sel, ]
    res$type <- sapply(Y, "[", 3L)[sel]

    i3 <- i2[sel]
    from <- i3 + 1L
    to <- c(i3[-1L] - 1L, n)

    x[i] <- gsub("^\t\t\t", "", x[i])
    Z <- mapply(get.product.others.gene, from, to)
    res$product <- Z[1L, ]
    gene <- Z[2L, ]
    res$others <- gsub("\t", ": ", Z[3L, ])
    if (!all(is.na(gene))) res$gene <- gene
    row.names(res) <- as.character(seq_len(nrow(res)))
    res
}

getAnnotationsGenBank <- function(access.nb, quiet = TRUE)
{
    N <- length(access.nb)
    res <- setNames(vector("list", N), access.nb)
    notfound <- NULL
    s1 <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id="
    s3 <- "&rettype=ft&retmode=text"
    for (i in 1:N) {
        if (!quiet) cat("\rDownloading annotations:", i, "/", N)
        URL <- paste0(s1, access.nb[i], s3)
        fl <- tempfile()
        ans <- try(download.file(URL, fl, quiet = TRUE), silent = TRUE)
        if (class(ans) == "try-error") {
            notfound <- c(notfound, access.nb[i])
            next
        }
        res[[i]] <- .parse.annotations.file(fl)
    }
    if (!quiet) cat(". Done.\n")
    if (!is.null(notfound)) {
        warning(paste0("cannot get features for the following accession(s):\n",
                       paste(notfound, collapse = ", ")))
        if (length(notfound) == N) return(NULL)
    }
    if (N == 1) res[[1L]] else res
}
