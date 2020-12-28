## apetools.R (2018-06-13)

##   APE Tools

## Copyright 2017-2018 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.file.extensions <-
    list(clustal = "aln", fasta = c("fas", "fasta", "fa"),
         fastq = c("fq", "fastq"), newick = c("nwk", "newick", "tre", "tree"),
         nexus = c("nex", "nexus"), phylip = "phy")

Xplorefiles <- function(from = "HOME", recursive = TRUE, ignore.case = TRUE)
{
    if (from == "HOME") from <- Sys.getenv("HOME")
    FILES <- list.files(path = from, recursive = recursive, full.names = TRUE)
    ext <- if (exists(".file.extensions", envir = .PlotPhyloEnv))
               get(".file.extensions", envir = .PlotPhyloEnv)
           else .file.extensions
    res <- vector("list", length(ext))
    names(res) <- names(ext)
    for (i in seq_along(res)) {
        e <- paste0("\\.", ext[[i]], "$")
        if (length(e) > 1) e <- paste(e, collapse = "|")
        x <- grep(e, FILES, ignore.case = ignore.case, value = TRUE)
        res[[i]] <- data.frame(File = x, Size = file.size(x),
                               stringsAsFactors = FALSE)
    }
    res
}

editFileExtensions <- function()
{
    foo <- function(x) {
        n <- length(x)
        if (n < m) x[(n + 1):m] <- NA
        x
    }
    res <- if (exists(".file.extensions", envir = .PlotPhyloEnv))
               get(".file.extensions", envir = .PlotPhyloEnv)
           else .file.extensions
    m <- max(lengths(res, FALSE))
    res <- lapply(res, foo)
    res <- as.data.frame(res, stringsAsFactors = FALSE)
    res <- edit(res)
    res <- lapply(res, function(x) x[!is.na(x)])
    assign(".file.extensions", res, envir = .PlotPhyloEnv)
}

bydir <- function(x)
{
    nofile <- which(sapply(x, nrow) == 0)
    if (length(nofile)) x <- x[-nofile]
    if (!length(x)) {
        cat("No file\n")
        return(invisible(NULL))
    }
    for (i in seq_along(x)) x[[i]]$Type <- names(x)[i]
    x <- do.call(rbind, x)
    x <- x[order(x$File), ]
    SPLIT <- strsplit(x$File, "/")
    LL <- lengths(SPLIT)
    foo <- function(i, PATH) {
        K <- grep(paste0("^", PATH, "/"), x$File)
        sel <- intersect(K, which(LL == i + 1L))
        if (length(sel)) {
            y <- x[sel, ]
            y$File <- gsub(".*/", "", y$File)
            cat("\n", PATH, "/\n", sep = "")
            print(y, row.names = FALSE)
        }
        if (length(sel) < length(K)) {
            d <- setdiff(K, sel)
            subdir <- unlist(lapply(SPLIT[d], "[", i + 1L))
            for (z in unique(subdir))
                foo(i + 1L, paste(PATH, z, sep = "/"))
        }
    }
    top <- unlist(lapply(SPLIT, "[", 1L))
    for (z in unique(top)) foo(1L, z)
}

Xplor <- function(from = "HOME")
{
    ext <- if (exists(".file.extensions", envir = .PlotPhyloEnv))
               get(".file.extensions", envir = .PlotPhyloEnv)
           else .file.extensions

    OUT <- paste0(tempfile(), ".html")
    mycat <- function(...) cat(..., sep = "", file = OUT, append = TRUE)

    FILES <- Xplorefiles(from = from)
    filetypes <- names(FILES)
    ## nb of files of each type:
    NR <- sapply(FILES, nrow)

    ## HTML header
    mycat('<html><title>Files Sorted by Type</title></head>')

    ## build the TOC
    mycat('<h2>File types searched:</h2>')
    mycat('<TABLE border=0>')
    mycat('<tr><th align="center"> Type </th><th align="center"> Number of files </th><th align="center"> Extensions* </th></tr>')
    for (type in filetypes) {
        mycat('<tr><td><a href=#', type, '>', type, '</a></td><td align="center">', NR[type],
                              '</td><td align="center">', paste(paste0(".", ext[[type]]), collapse = " "), '</td></tr>')
    }
    mycat('</TABLE>')

    mycat('<br>*Case-independent<br>To change the files extensions, type in R: <font face = "courier">editFileExtensions()</font><br>')

    if (all(NR == 0)) {
        browseURL(OUT)
        return(invisible(NULL))
    }

    OUTBYDIR <- paste0(tempfile(), ".html")
    sink(OUTBYDIR)
    cat('<html><title>Files Sorted by Directory</title></head>')
    .bydir.html(FILES)
    cat('</html>')
    sink(NULL)
    mycat('<br><a target="blank" href=', OUTBYDIR, '><h2>Files sorted by directory (in new tab)</h2></a><br>')

    for (type in filetypes) {
        nr <- NR[type]
        mycat('<a name=', type, '></a><h1>', toupper(type), '</h1>')
        if (nr == 0) {
            mycat('no file of this type')
            next
        }
        DF <- FILES[[type]]
        mycat('<TABLE border=1>')
        mycat('<tr><th> File name </th><th align="center"> Size (KB) </th></tr>')
        for (i in 1:nr)
            mycat('<tr><td><a href="file://', DF[i, 1], '">', DF[i, 1],
                  '</a></td><td align="right">', round(DF[i, 2]/1000, 1), '</td></tr>')
        mycat('</TABLE>')
    }
    mycat('</html>')
    browseURL(OUT)
}

.bydir.html <- function(x)
{
    nofile <- which(sapply(x, nrow) == 0)
    if (length(nofile)) x <- x[-nofile]
    if (!length(x)) return(NULL)
    for (i in seq_along(x)) x[[i]]$Type <- names(x)[i]
    x <- do.call(rbind, x)
    x <- x[order(x$File), ]
    SPLIT <- strsplit(x$File, "/")
    LL <- lengths(SPLIT)
    foo <- function(i, PATH) {
        K <- grep(paste0("^", PATH, "/"), x$File)
        sel <- intersect(K, which(LL == i + 1L))
        if (length(sel)) {
            y <- x[sel, ]
            y$File <- gsub(".*/", "", y$File)
            cat('<h2><a href="file://', PATH, '">', PATH, '/</a></h2>', sep = "")
            cat('<TABLE border=1>')
            cat('<tr><th> File </th><th align="center"> Size (KB) </th><th align="center"> Type </th></tr>')
            for (i in 1:nrow(y))
                cat('<tr><td><a href="file://', PATH, '/', y[i, 1], '">', y[i, 1],
                    '</a></td><td align="right">', round(y[i, 2]/1000, 1),
                    '</td><td align="right">', y[i, 3], '</td></tr>', sep = "")
            cat('</TABLE><br>')
        }
        if (length(sel) < length(K)) {
            d <- setdiff(K, sel)
            subdir <- unlist(lapply(SPLIT[d], "[", i + 1L))
            for (z in unique(subdir))
                foo(i + 1L, paste(PATH, z, sep = "/"))
        }
    }
    top <- unlist(lapply(SPLIT, "[", 1L))
    for (z in unique(top)) foo(1L, z)
}
