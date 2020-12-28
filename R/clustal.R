## clustal.R (2018-03-28)

##   Multiple Sequence Alignment with External Applications

## Copyright 2011-2018 Emmanuel Paradis, 2018 Franz Krah

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.errorAlignment <- function(exec, prog)
{
    dirs <- strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
    paste0("\n   cannot find executable ", sQuote(exec), " on your computer.\n",
           "  It is recommended that you place the executable of ", prog, "\n",
           "  in a directory on the PATH of your computer which is:\n",
           paste(sort(dirs), collapse = "\n"))
}

clustalomega <- function (x, y, guide.tree, exec = NULL, MoreArgs = "",
                           quiet = TRUE, original.ordering = TRUE, file)
{
    os <- Sys.info()[1]
    if (is.null(exec)) {
        exec <- switch(os, Linux = "clustalo", Darwin = "clustalo",
                       Windows = "clustalo.exe")
    }
    if (missing(x)) {
        out <- system(paste(exec, "-h"))
        if (out == 127)
            stop(.errorAlignment(exec, "Clustal-Omega"))
        return(invisible(NULL))
    }
    type <- if (inherits(x, "DNAbin")) "DNA" else "AA"
    if (type == "AA" && !inherits(x, "AAbin"))
        stop("'x' should be of class \"DNAbin\" or \"AAbin\"")
    noy <- missing(y)

    fns <- character(4)
    for (i in 1:3)
        fns[i] <- tempfile(pattern = "clustal", tmpdir = tempdir(), fileext = ".fas")
    fns[4] <-  tempfile(pattern = "guidetree", tmpdir = tempdir(), fileext = ".nwk")
    unlink(fns[file.exists(fns)])

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))
    write.FASTA(x, fns[1])

    if (noy) {
        opts <- paste("-i", fns[1], "-o", fns[3], "--force")

        ## add input guide tree
        if (!missing(guide.tree)) {
            if (!inherits(guide.tree, "phylo"))
                stop("object 'guide.tree' is not of class \"phylo\"")
            if (length(setdiff(labels.bak, guide.tree$tip.label)))
                stop("guide tree does not match sequence names")
            guide.tree$tip.label[match(guide.tree$tip.label, labels.bak)] <- names(x)
            if (!is.binary(guide.tree)) guide.tree <- multi2di(guide.tree)
            if (is.null(guide.tree$edge.length))
                guide.tree$edge.length <- rep(1, Nedge(guide.tree))
            write.tree(guide.tree, fns[4])
            opts <- paste(opts, paste("--guidetree-in", fns[4]))
        }
    } else {
        y <- as.list(y)
        labels.baky <- names(y)
        names(y) <- paste0("Id", length(x) + 1:length(y))
        write.FASTA(y, fns[2])
        if (length(y) == 1) {
            opts <- paste("-i", fns[1],"--profile1",
                          fns[2], "-o", fns[3], "--force")
        } else {
            opts <- paste("--profile1", fns[1],"--profile2",
                          fns[2], "-o", fns[3], "--force")
        }
    }

    opts <- paste(opts, MoreArgs)
    if (!quiet) opts <- paste(opts, "-v")
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
        stop(.errorAlignment(exec, "Clustal-Omega"))
    res <- as.matrix(read.FASTA(fns[3], type))

    if (noy) {
        if (original.ordering)
            res <- res[labels(x), ]
        rownames(res) <- labels.bak
    } else {
        if (original.ordering)
            res <- res[c(labels(x), labels(y)), ]
        rownames(res) <- c(labels.bak, labels.baky)
    }

    unlink(fns[file.exists(fns)])
    if (missing(file)) return(res) else write.FASTA(res, file)
}

clustal <-
    function(x, y, guide.tree, pw.gapopen = 10, pw.gapext = 0.1, gapopen = 10,
             gapext = 0.2, exec = NULL, MoreArgs = "", quiet = TRUE,
             original.ordering = TRUE, file)
{
    os <- Sys.info()[1]
    if (is.null(exec)) {
        exec <- switch(os, Linux = "clustalw", Darwin = "clustalw2",
                       Windows = "clustalw2.exe")
    }
    if (missing(x)) {
        out <- system(paste(exec, "-help"))
        if (out == 127)
            stop(.errorAlignment(exec, "Clustal"))
        return(invisible(NULL))
    }
    type <- if (inherits(x, "DNAbin")) "DNA" else "AA"
    if (type == "AA" && !inherits(x, "AAbin"))
        stop("'x' should be of class \"DNAbin\" or \"AAbin\"")
    noy <- missing(y)

    fns <- character(4)
    for (i in 1:3)
        fns[i] <- tempfile(pattern = "clustal", tmpdir = tempdir(), fileext = ".fas")
    fns[4] <-  tempfile(pattern = "guidetree", tmpdir = tempdir(), fileext = ".nwk")
    unlink(fns[file.exists(fns)])

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))
    write.FASTA(x, fns[1])

    if (noy) {
        prefix <- c("-INFILE", "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN","-GAPEXT", "-OUTFILE","-OUTPUT")
        suffix <- c(fns[1], pw.gapopen, pw.gapext, gapopen, gapext, fns[3], "FASTA")

        ## add input guide tree
        if (!missing(guide.tree)) {
            if (!inherits(guide.tree, "phylo"))
                stop("object 'guide.tree' is not of class \"phylo\"")
            if (length(setdiff(labels.bak, guide.tree$tip.label)))
                stop("guide tree does not match sequence names")
            guide.tree$tip.label[match(guide.tree$tip.label, labels.bak)] <- names(x)
            if (!is.binary(guide.tree)) guide.tree <- multi2di(guide.tree)
            if (is.null(guide.tree$edge.length))
                guide.tree$edge.length <- rep(1, Nedge(guide.tree))
            write.tree(guide.tree, fns[4])
            prefix <- c(prefix, "-USETREE")
            suffix <- c(suffix, fns[4])
        }
    } else {
        y <- as.list(y)
        labels.baky <- names(y)
        names(y) <- paste0("Id", length(x) + 1:length(y))
        write.FASTA(y, fns[2])
        prefix <- c("-PROFILE1", "-PROFILE2", "-PWGAPOPEN", "-PWGAPEXT",
                    "-GAPOPEN","-GAPEXT", "-OUTFILE","-OUTPUT")
        suffix <- c(fns[1], fns[2], pw.gapopen, pw.gapext,
                    gapopen, gapext, fns[3], "FASTA")

    }

    opts <- paste(prefix, suffix, sep = "=", collapse = " ")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
        stop(.errorAlignment(exec, "Clustal"))
    res <- as.matrix(read.FASTA(fns[3], type))

    if (noy) {
        if (original.ordering)
            res <- res[labels(x), ]
        rownames(res) <- labels.bak
    } else {
        if (original.ordering)
            res <- res[c(labels(x), labels(y)), ]
        rownames(res) <- c(labels.bak, labels.baky)
    }

    unlink(fns[file.exists(fns)])
    if (missing(file)) return(res) else write.FASTA(res, file)
}

muscle <- function (x, y, guide.tree, exec = "muscle", MoreArgs = "",
                     quiet = TRUE, original.ordering = TRUE, file)
{
    if (missing(x)) {
        out <- system(exec)
        if (out == 127) stop(.errorAlignment(exec, "MUSCLE"))
        return(invisible(NULL))
    }
    type <- if (inherits(x, "DNAbin")) "DNA" else "AA"
    if (type == "AA" && !inherits(x, "AAbin"))
        stop("'x' should be of class \"DNAbin\" or \"AAbin\"")
    noy <- missing(y)

    ## Produce TEMP files
    fns <- character(4)
    for (i in 1:3)
        fns[i] <- tempfile(pattern = "muscle", tmpdir = tempdir(), fileext = ".fas")
    fns[4] <- tempfile(pattern = "guidetree", tmpdir = tempdir(), fileext = ".nwk")
    unlink(fns[file.exists(fns)])

    ## Write input sequences x to file
    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))
    write.FASTA(x, fns[1])

    ## Run muscle for X
    if (noy) {
        opts <- paste("-in", fns[1], "-out", fns[3])
        ## add input guide tree
        if (!missing(guide.tree)) {
            if (!inherits(guide.tree, "phylo"))
                stop("object 'guide.tree' is not of class \"phylo\"")
            if (length(setdiff(labels.bak, guide.tree$tip.label)))
                stop("guide tree does not match sequence names")
            guide.tree$tip.label[match(guide.tree$tip.label, labels.bak)] <- names(x)
            if (!is.binary(guide.tree)) guide.tree <- multi2di(guide.tree)
            if (is.null(guide.tree$edge.length))
                guide.tree$edge.length <- rep(1, Nedge(guide.tree))
            write.tree(guide.tree, fns[4])
            opts <- paste(opts, paste("-usetree_nowarn", fns[4]))
        }
    } else {
        y <- as.list(y)
        labels.baky <- names(y)
        names(y) <- paste0("Id", length(x) + 1:length(y))
        write.FASTA(y, fns[2])
        opts <- paste("-profile", "-in1", fns[1],"-in2", fns[2], "-out", fns[3])
    }

    if (quiet) opts <- paste(opts, "-quiet")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts))
    if (out == 127)
        stop(.errorAlignment(exec, "MUSCLE"))
    res <- as.matrix(read.FASTA(fns[3], type))

    if (noy) {
        if (original.ordering)
            res <- res[labels(x), ]
        rownames(res) <- labels.bak
    } else {
        if (original.ordering)
            res <- res[c(labels(x), labels(y)), ]
        rownames(res) <- c(labels.bak, labels.baky)
    }

    unlink(fns[file.exists(fns)])
    if (missing(file)) return(res) else write.FASTA(res, file)
}

tcoffee <- function(x, exec = "t_coffee", MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
    if (missing(x)) {
        out <- system(exec)
        if (out == 127) stop(.errorAlignment(exec, "T-Coffee"))
        return(invisible(NULL))
    }

    x <- as.list(x)
    labels.bak <- names(x)
    names(x) <- paste0("Id", 1:length(x))

    d <- tempdir()
    od <- setwd(d)
    on.exit(setwd(od))
    inf <- "input_tcoffee.fas"
    write.dna(x, inf, "fasta")
    opts <- paste(inf, MoreArgs)
    if (quiet) opts <- paste(opts, "-quiet=nothing")
    out <- system(paste(exec, opts))
    if (out == 127) stop(.errorAlignment(exec, "T-Coffee"))
    res <- read.dna("input_tcoffee.aln", "clustal")
    if (original.ordering) res <- res[labels(x), ]
    rownames(res) <- labels.bak
    res
}
