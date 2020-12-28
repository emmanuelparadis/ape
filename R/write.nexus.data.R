## write.nexus.data.R (2018-06-23)

##   Write Character Data in NEXUS Format

## Copyright 2006-2015 Johan Nylander, Emmanuel Paradis, 2018 Thomas Guillerme

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

write.nexus.data <-
    function(x, file, format = "dna", datablock = TRUE,
             interleaved = TRUE, charsperline = NULL,
             gap = NULL, missing = NULL)
{
### TODO: Standard data, mixed data, nice indent

    format <- match.arg(toupper(format), c("DNA", "PROTEIN", "STANDARD", "CONTINUOUS"))

    if (inherits(x, "DNAbin") && format != "DNA") {
        format <- "DNA"
        warning("object 'x' is of class DNAbin: format forced to DNA")
    }
    if (inherits(x, "AAbin") && format != "PROTEIN") {
        format <- "PROTEIN"
        warning("object 'x' is of class AAbin: format forced to PROTEIN")
    }

    indent          <- "  "  # Two blanks
    maxtax          <- 5     # Max nr of taxon names to be printed on a line
    defcharsperline <- 80    # Default nr of characters per line if interleaved
    defgap          <- "-"   # Default gap character
    defmissing      <- "?"   # Default missing data character

    if (is.matrix(x)) {
        if (inherits(x, "DNAbin")) x <- as.list(x) else {
            xbak <- x
            x <- vector("list", nrow(xbak))
            for (i in seq_along(x)) x[[i]] <- xbak[i, ]
            names(x) <- rownames(xbak)
            rm(xbak)
        }
    }

    ntax <- length(x)
    nchars <- length(x[[1]])

    zz <- file(file, "w")

    if (is.null(names(x))) names(x) <- as.character(1:ntax)

    fcat <- function(..., file = zz)
        cat(..., file = file, sep = "", append = TRUE)

    find.max.length <- function(x) max(nchar(x))

    print.matrix <- function(x, dindent = "    ", collapse = "") {
        Names <- names(x)
        printlength <- find.max.length(Names) + 2
        if (!interleaved) {
            for (i in seq_along(x)) {
                sequence <- paste(x[[i]], collapse = collapse)
                taxon <- Names[i]
                thestring <- sprintf("%-*s%s%s", printlength, taxon, dindent, sequence)
                fcat(indent, indent, thestring, "\n")
            }
        } else {
            ntimes <- ceiling(nchars/charsperline)
            start <- 1
            end <- charsperline
            for (j in seq_len(ntimes)) {
                for (i in seq_along(x)) {
                    sequence <- paste(x[[i]][start:end], collapse = collapse)
                    taxon <- Names[i]
                    thestring <- sprintf("%-*s%s%s", printlength, taxon, dindent, sequence)
                    fcat(indent, indent, thestring, "\n")
                }
                if (j < ntimes) fcat("\n")
                start <- start + charsperline
                end <- end + charsperline
                if (end > nchars) end <- nchars
            }
        }
    }

    if (inherits(x, "DNAbin") || inherits(x, "AAbin")) x <- as.character(x)

    fcat("#NEXUS\n[Data written by write.nexus.data.R, ", date(), "]\n")

    NCHAR <- paste("NCHAR=", nchars, sep = "")
    NTAX <- paste0("NTAX=", ntax)
    DATATYPE <- paste0("DATATYPE=", format) # fix by Robin Cristofari (2015-02-04)

    if (is.null(charsperline)) {
        if (nchars <= defcharsperline) {
            charsperline <- nchars
            interleaved <- FALSE
        } else charsperline <- defcharsperline
    }

    if (is.null(missing)) missing <- defmissing
    MISSING <- paste0("MISSING=", missing)

    if (is.null(gap)) gap <- defgap
    GAP <- paste0("GAP=", gap)

    INTERLEAVE <- if (interleaved) "INTERLEAVE=YES" else "INTERLEAVE=NO"

    if (datablock) {
        fcat("BEGIN DATA;\n")
        fcat(indent, "DIMENSIONS ", NTAX, " ", NCHAR, ";\n")

        ## <FIXME> only DNA and PROTEIN is supported for the moment, so the
        ## following 'if' is not needed
        ## if (format %in% c("DNA", "PROTEIN")) # from Francois Michonneau (2009-10-02)
        if(format != "STANDARD") {
            fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", GAP, " ", INTERLEAVE, ";\n")
        } else {
            fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", GAP, " ", INTERLEAVE, " symbols=\"0123456789\";\n")
        }
        ## </FIXME>

        fcat(indent, "MATRIX\n")
        if(format != "CONTINUOUS") {
            print.matrix(x)
        } else {
            print.matrix(x, collapse = "\t")
        }
        fcat(indent, ";\nEND;\n\n")
    } else {
        fcat("BEGIN TAXA;\n")
        fcat(indent, "DIMENSIONS", " ", NTAX, ";\n")
        fcat(indent, "TAXLABELS\n")
        fcat(indent, indent)
        j <- 0
        for (i in seq_len(ntax)) {
            fcat(names(x[i]), " ")
            j <- j + 1
            if (j == maxtax) {
                fcat("\n", indent, indent)
                j <- 0
            }
        }
        fcat("\n", indent, ";\n")
        fcat("END;\n\nBEGIN CHARACTERS;\n")
        fcat(indent, "DIMENSIONS", " ", NCHAR, ";\n")

        ## <FIXME> only DNA and PROTEIN is supported for the moment, so the
        ## following 'if' is not needed
        ## if (format %in% c("DNA", "PROTEIN"))
        if(format != "STANDARD") {
            fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", DATATYPE, " ", INTERLEAVE, ";\n")
        } else {
            fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", DATATYPE, " ", INTERLEAVE, " symbols=\"0123456789\";\n")
        }
        ## </FIXME>

        fcat(indent,"MATRIX\n")
        if(format != "CONTINUOUS") {
            print.matrix(x)
        } else {
            print.matrix(x, collapse = "\t")
        }
        fcat(indent, ";\nEND;\n\n")
    }
    close(zz)
}
