## write.nexus.R (2017-09-08)

##   Write Tree File in Nexus Format

## Copyright 2003-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

write.nexus <- function(..., file = "", translate = TRUE)
{
    obj <- .getTreesFromDotdotdot(...)
    ntree <- length(obj)
    cat("#NEXUS\n", file = file)
    cat(paste("[R-package APE, ", date(), "]\n\n", sep = ""),
        file = file, append = TRUE)

    N <- length(obj[[1]]$tip.label)

    cat("BEGIN TAXA;\n", file = file, append = TRUE)
    cat(paste("\tDIMENSIONS NTAX = ", N, ";\n", sep = ""),
        file = file, append = TRUE)
    cat("\tTAXLABELS\n", file = file, append = TRUE)
    cat(paste("\t\t", obj[[1]]$tip.label, sep = ""),
        sep = "\n", file = file, append = TRUE)
    cat("\t;\n", file = file, append = TRUE)
    cat("END;\n", file = file, append = TRUE)

    cat("BEGIN TREES;\n", file = file, append = TRUE)
    if (translate) {
        cat("\tTRANSLATE\n", file = file, append = TRUE)
        obj <- .compressTipLabel(obj)
        X <- paste("\t\t", 1:N, "\t", attr(obj, "TipLabel"), ",", sep = "")
        ## We remove the last comma:
        X[length(X)] <- gsub(",", "", X[length(X)])
        cat(X, file = file, append = TRUE, sep = "\n")
        cat("\t;\n", file = file, append = TRUE)
        class(obj) <- NULL
        for (i in 1:ntree)
            obj[[i]]$tip.label <- as.character(1:N)
    } else {
        if (is.null(attr(obj, "TipLabel"))) {
            for (i in 1:ntree)
                obj[[i]]$tip.label <- checkLabel(obj[[i]]$tip.label)
        } else {
            attr(obj, "TipLabel") <- checkLabel(attr(obj, "TipLabel"))
            obj <- .uncompressTipLabel(obj)
        }
    }

    title <- names(obj)
    if (is.null(title))
        title <- rep("UNTITLED", ntree)
    else {
        if (any(s <- title == "")) title[s] <- "UNTITLED"
    }

    for (i in 1:ntree) {
        if (class(obj[[i]]) != "phylo") next
        root.tag <- if (is.rooted(obj[[i]])) "= [&R] " else "= [&U] "
        cat("\tTREE *", title[i], root.tag, file = file, append = TRUE)
        cat(write.tree(obj[[i]], file = ""),
            "\n", sep = "", file = file, append = TRUE)
    }
    cat("END;\n", file = file, append = TRUE)
}
