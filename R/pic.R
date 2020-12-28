## pic.R (2017-08-22)

##   Phylogenetically Independent Contrasts

## Copyright 2002-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

pic <- function(x, phy, scaled = TRUE, var.contrasts = FALSE, rescaled.tree = FALSE)
{
    if (!inherits(phy, "phylo")) stop("object 'phy' is not of class \"phylo\"")
    if (is.null(phy$edge.length)) stop("your tree has no branch lengths")
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if (nb.node != nb.tip - 1)
        stop("'phy' is not rooted and fully dichotomous")
    if (length(x) != nb.tip)
        stop("length of phenotypic and of phylogenetic data do not match")
    if (any(is.na(x)))
        stop("missing data in 'x': you may consider removing the species with missing data from your tree with the function 'drop.tip'.")

    phy <- reorder(phy, "postorder")
    phenotype <- numeric(nb.tip + nb.node)

    if (is.null(names(x))) {
        phenotype[1:nb.tip] <- x
    } else {
        if (all(names(x) %in% phy$tip.label))
          phenotype[1:nb.tip] <- x[phy$tip.label]
        else {
            phenotype[1:nb.tip] <- x
            warning("the names of argument 'x' and the tip labels of the tree did not match: the former were ignored in the analysis.")
        }
    }

    ## No need to copy the branch lengths: they are rescaled
    ## in the C code, so it's important to leave the default
    ## `DUP = TRUE' of .C.
    ans <- .C(C_pic, as.integer(nb.tip),
              as.integer(phy$edge[, 1]), as.integer(phy$edge[, 2]),
              as.double(phy$edge.length), as.double(phenotype),
              double(nb.node), double(nb.node),
              as.integer(var.contrasts), as.integer(scaled))

    contr <- ans[[6]]
    lbls <-
        if (is.null(phy$node.label)) as.character(1:nb.node + nb.tip)
        else phy$node.label
    if (var.contrasts) {
        contr <- cbind(contr, ans[[7]])
        dimnames(contr) <- list(lbls, c("contrasts", "variance"))
    } else names(contr) <- lbls
    if (rescaled.tree) {
        phy$edge.length <- ans[[4]]
        contr <- list(contr = contr, rescaled.tree = phy)
    }
    contr
}

pic.ortho <- function(x, phy, var.contrasts = FALSE, intra = FALSE)
{
    n <- length(x)
    m <- n - 1L # number of nodes
    phy <- reorder(phy, "postorder")
    xx <- unlist(lapply(x, mean)) # 'x' in Felsenstein's paper
    xx <- c(xx, numeric(m))
    delta.v <- numeric(n + m)
    s <- 1/lengths(x)
    s <- c(s, numeric(m))
    contrast <- var.cont <- numeric(m)

    i <- 1L
    while (i < m + n) {
        d1 <- phy$edge[i, 2]
        d2 <- phy$edge[i + 1L, 2]
        a <- phy$edge[i, 1]
        tmp1 <- 1/(phy$edge.length[i] + delta.v[d1])
        tmp2 <- 1/(phy$edge.length[i + 1L] + delta.v[d2])
        xx[a] <- (tmp1 * xx[d1] + tmp2 * xx[d2])/(tmp1 + tmp2)
        delta.v[a] <- 1/(tmp1 + tmp2)
        f1 <- tmp1/(tmp1 + tmp2)
        f2 <- tmp2/(tmp1 + tmp2)
        s[a] <- f1*f1 * s[d1] + f2*f2 * s[d2]
        tmp <- 1/(s[d1] + s[d2])
        contrast[a - n] <- (xx[d1] - xx[d2]) * sqrt(tmp)
        var.cont[a - n] <- (1/tmp1 + 1/tmp2) * tmp
        i <- i + 2L
    }

    lbls <-
        if (is.null(phy$node.label)) as.character(1:m + n)
        else phy$node.label

    if (var.contrasts) {
        contrast <- cbind(contrast, var.cont)
        dimnames(contrast) <- list(lbls, c("contrasts", "variance"))
    } else names(contrast) <- lbls

    if (intra) {
        intraspe.ctr <- function(x) {
            k <- length(x) - 1L
            if (!k) return(NULL)
            ctr <- numeric(k)
            ctr[1L] <- x[1L] - x[2L]
            if (k > 1)
                for (i in 2:k)
                    ctr[i] <- x[i + 1L] - mean(x[1:i])
            sqrt((1:k)/(1:k + 1)) * ctr
        }
        tmp <- lapply(x, intraspe.ctr)
        names(tmp) <- phy$tip.label
        attr(contrast, "intra") <- tmp
    }

    contrast
}

varCompPhylip <- function(x, phy, exec = NULL)
{
    n <- Ntip(phy)
    if (is.vector(x)) x <- as.list(x)
    if (is.matrix(x) || is.data.frame(x)) {
        tmpx <- vector("list", n)
        for (i in 1:n) tmpx[[i]] <- x[i, , drop = FALSE]
        names(tmpx) <- rownames(x)
        x <- tmpx
    }
    p <- if (is.vector(x[[1]])) 1L else ncol(x[[1]])
    if (!is.null(names(x))) x <- x[phy$tip.label]

    phy <- makeLabel(phy, len = 10)
    lbs <- phy$tip.label

    ni <- sapply(x, function(xx) if (is.vector(xx)) 1L else nrow(xx))

    pfx <- tempdir()
    write.tree(phy, file = paste(pfx, "intree", sep = "/"))
    infile <- paste(pfx, "infile", sep = "/")
    file.create(infile)
    cat(n, " ", p, "\n", sep = "", file = infile, append = TRUE)
    for (i in 1:n) {
        cat(lbs[i], file = infile, append = TRUE)
        ## can surely be better but OK for the moment:
        cat(paste(rep(" ", 11 - nchar(lbs[i])), collapse = ""),
            file = infile, append = TRUE)
        cat(ni[i], "\n", sep = "", file = infile, append = TRUE)
        if (ni[i] == 1) {
            cat(x[[i]], sep = " ", file = infile, append = TRUE)
            cat("\n", file = infile, append = TRUE)
        } else write(t(x[[i]]), file = infile, ncolumns = p, append = TRUE)
    }

    if (is.null(exec))
        exec <-
            if (.Platform$OS.type == "unix") "phylip contrast"
            else "contrast"

    odir <- setwd(pfx)
    on.exit(setwd(odir))
    if (file.exists("outfile")) unlink("outfile")
    system(exec, intern = TRUE, input = c("W", "A", "Y"))
    varA <- scan("outfile", skip = 7, nlines = p, quiet = TRUE)
    varE <- scan("outfile", skip = 11 + p, nlines = p, quiet = TRUE)
    if (p > 1) {
        varA <- matrix(varA, p, p, byrow = TRUE)
        varE <- matrix(varE, p, p, byrow = TRUE)
    }
    list(varA = varA, varE = varE)
}
