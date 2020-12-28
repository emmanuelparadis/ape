## read.gff.R (2018-05-22)

##   Read GFF Files

## Copyright 2016-2018 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

read.gff <- function(file, na.strings = c(".", "?"), GFF3 = TRUE)
{
    w <- list("", "", "", 0L, 0L, 0, "", "", "")
    x <- scan(file, w, sep = "\t", quote = "", quiet = TRUE,
              na.strings = na.strings, comment.char = "#")
    for (i in c(1, 2, 3, 7, 8)) x[[i]] <- factor(x[[i]])
    names(x) <- c("seqid", "source", "type", "start", "end",
                  "score", "strand", "phase", "attributes")
    if (!GFF3) {
        names(x) <- c("seqname", "source", "feature", "start", "end",
                      "score", "strand", "frame", "attributes")
    }
    n <- length(x[[1]])
    attr(x, "row.names") <- as.character(seq_len(n))
    class(x) <- "data.frame"
    x
}
