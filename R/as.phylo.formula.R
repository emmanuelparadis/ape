## as.phylo.formula.R (2018-09-17)

##   Conversion from Taxonomy Variables to Phylogenetic Trees

## Copyright 2005-2018 Julien Dutheil, 2018 Eric Marcon

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

as.phylo.formula <- function(x, data = parent.frame(), collapse = TRUE, ...)
{
    ## Testing formula syntax:
    err <- "Formula must be of the kind ~A1/A2/.../An."
    if (length(x) != 2) stop(err)
    if (x[[1]] != "~") stop(err)
    f <- x[[2]]
    taxo <- list()
    while (length(f) == 3) {
        if (f[[1]] != "/") stop(err)
        f3.txt <- deparse(f[[3]])
        if (!is.factor(data[[f3.txt]]))
            stop(paste("Variable", f3.txt, "must be a factor"))
        taxo[[f3.txt]] <- data[[f3.txt]]
        if (length(f) > 1) f <- f[[2]]
    }
    f.txt <- deparse(f)
    if (!is.factor(data[[f.txt]]))
        stop(paste("Variable", f.txt, "must be a factor."))
    taxo[[f.txt]] <- data[[f.txt]]
    taxo.data <- as.data.frame(taxo)
    leaves.names <- as.character(taxo.data[, 1])
    taxo.data[, 1] <- 1:nrow(taxo.data)
    ## Now builds the phylogeny:
    f.rec <- function(subtaxo) { # Recurrent utility function
        u <- ncol(subtaxo)
        levels <- unique(subtaxo[,u])
        if (u == 1) {
            if (length(levels) != nrow(subtaxo))
                warning("leaves names are not unique.")
            return(as.character(subtaxo[, 1]))
        }
        t <- character(length(levels))
        for (l in 1:length(levels)) {
            x <- f.rec(subtaxo[subtaxo[,u] == levels[l], ][1:(u - 1)])
            t[l] <- paste0("(", paste(x, collapse=","), ")")
        }
        t
    }
    string <- paste0("(", paste(f.rec(taxo.data), collapse = ","), ");")
    phy <- read.tree(text = string)
    if (collapse) phy <- collapse.singles(phy)
    phy$tip.label <- leaves.names[as.numeric(phy$tip.label)]
    phy
}
