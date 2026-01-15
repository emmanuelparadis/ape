## read.nexus.R (2026-01-14)

##   Read Tree File in Nexus Format

## Copyright 2003-2025 Emmanuel Paradis and 2010-2025 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.treeBuild <- function(x)
{
    if (!length(grep(",", x))) {
        phy <- list(edge = matrix(c(2L, 1L), 1, 2), Nnode = 1L)
        x <- unlist(strsplit(x, "[\\(\\):;]"))
        phy$tip.label <- x[2]
        phy$edge.length <- as.numeric(x[3])
        phy$node.label <- x[4]
    } else {
        phy <- .Call(treeBuild, x)
        dim(phy[[1]]) <- c(length(phy[[1]])/2, 2)
        nms <- c("edge", "edge.length", "Nnode", "node.label", "tip.label", "root.edge")
        if (length(phy) == 5) nms <- nms[-6]
        names(phy) <- nms
    }
    if (all(phy$node.label == "")) phy$node.label <- NULL
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

.cladoBuild <- function(x)
{
    if (!length(grep(",", x))) {
        ## only one tip but can be several nodes (GH's issue #104)
        Nnode <- length(gregexpr("\\)", x)[[1]])
        edge <- if (Nnode == 1L) 2:1 else c(2L, rep(3:(Nnode + 1L), each = 2), 1L)
        dim(edge) <- c(Nnode, 2L)
        phy <- list(edge = edge, Nnode = Nnode)
        labs <- unlist(strsplit(x, "[\\(\\);]"))
        phy$tip.label <- labs[Nnode + 1L]
        phy$node.label <- labs[(Nnode + 2L):length(labs)]
    } else {
        phy <- .Call(cladoBuild, x)
        dim(phy[[1]]) <- c(length(phy[[1]])/2, 2)
        nms <- c("edge", "Nnode", "node.label", "tip.label", "root.edge")
        if (length(phy) == 4) nms <- nms[-5]
        names(phy) <- nms
    }
    if (all(phy$node.label == "")) phy$node.label <- NULL
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

.evonetBuild <- function(x){ # , colon=TRUE){
    colon <- grepl(":", x)
    if (colon) {
        info <- extract_hybrid_info(x)
        for (i in seq_along(info$id)) x <- sub(info$name[i], info$id[i], x)
        y <- .treeBuild(x)
        z <- as.evonet(y, info = info)
    } else {
        y <- .cladoBuild(x)
        z <- as.evonet(y)
    }
    z
}

.decodeTRANSLATE <- function(x)
{
    z <- paste(x, collapse = " ")
    ## delete the ending semicolon:
    z <- sub("[[:space:]]*;[[:space:]]*$", "", z)
    z <- gsub("^[[:space:]]*TRANSLATE[[:space:]]*", "",
              z, ignore.case = TRUE)
    dq <- gregexpr('"', z)[[1]] # find double quotes
    sq <- gregexpr("'", z)[[1]] # find single quotes
    DQ <- SQ <- FALSE
    if (dq[1] > -1) {
        DQ <- TRUE
        if (!length(dq) %% 2)
            warning("found an odd number of double quotes in TRANSLATE command of NEXUS file")
    }
    if (sq[1] > -1) {
        SQ <- TRUE
        if (!length(sq) %% 2)
            warning("found an odd number of single quotes in TRANSLATE command of NEXUS file")
    }
    s <- c(TRUE, FALSE)
    if (DQ && SQ) {
        nDQ <- length(DQ)
        nSQ <- length(SQ)
        if (DQ[1] < SQ[1] && DQ[nDQ] > SQ[nSQ]) {
            start <- DQ[s]; stop <- DQ[!s]
            SQ <- FALSE
        } else {
            if (DQ[1] > SQ[1] && DQ[nDQ] < SQ[nSQ]) {
                start <- SQ[s]; stop <- SQ[!s]
                DQ <- FALSE
            }
        }
    }
    if (DQ && SQ)
        stop("inconsistent use of single and double quotes in TRANSLATE command of NEXUS file")
    if (DQ || SQ) {
        Z <- list(x = z)
        n <- length(stop)
        labs <- mapply(substr, start, stop, MoreArgs = Z)
        start2 <- c(1L, stop[-n] + 1L)
        stop2 <- start[-n] - 1L
        tokens <- mapply(substr, start2, stop2, MoreArgs = Z)
        tokens <- gsub("[[:space:]]*,*[[:space:]]*", "", tokens)
        if (grepl("[:space:]", tokens))
            warning("spaces found in tokens of NEXUS file")
        res <- cbind(tokens, labs)
    } else { # if no quotes at all
        res <- unlist(strsplit(z, "[[:space:]]*,[[:space:]]*"))
        res <- unlist(strsplit(res, "[[:space:]]+"))
        res <- matrix(res, ncol = 2, byrow = TRUE)
    }
    if (anyDuplicated(res))
        warning("found duplicates among tokens and labels")
    res
}

.translateTOKENS <- function(x, TABLE)
{
    ## 1) check the tokens:
    m <- match(x, TABLE[, 1L])
    NAs <- is.na(m)
    if (!any(NAs)) return(TABLE[m, 2L])
    i <- which(!NAs)
    if (length(i)) # in case no match at all
        x[i] <- TABLE[m[i], 2L]
    ## 2) check the labels:
    m <- match(x, TABLE[, 2L])
    i <- which(is.na(m))
    if (!length(i)) return(x)
    ## 3) check if this can be a taxon number....
    suppressWarnings(test <- as.integer(x[i]))
    valid <- !is.na(test) & test <= length(x)
    if (any(valid))
        x[i][valid] <- TABLE[test[valid], 2L]
    x
}

read.nexus <- function(file, tree.names = NULL, force.multi = FALSE)
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    ## remove all comments
    ## (this might not work if there are square brackets within the comments)
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) { # in case there are no comments at all this block is skipped
        w <- LEFT == RIGHT
        if (any(w)) { # in case all comments use at least 2 lines
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
            ## The above regexp was quite tough to find: it makes
            ## possible to delete series of comments on the same line:
            ##       ...[...]xxx[...]...
            ## without deleting the "xxx". This regexp is in three parts:
            ##       \\[   [^]]*   \\]
            ## where [^]]* means "any character, except "]", repeated zero
            ## or more times" (note that the ']' is not escaped here).
            ## The previous version was:
            ##       X[s] <- gsub("\\[.*\\]", "", X[s])
            ## which deleted the "xxx". (EP  2008-06-24)
        }
        w <- !w
        if (any(w)) {
            s <- LEFT[w]
            X[s] <- gsub("\\[.*", "", X[s])
            sb <- RIGHT[w]
            X[sb] <- gsub(".*\\]", "", X[sb])
            if (any(s < sb - 1))
                X <- X[-unlist(mapply(":", (s + 1), (sb - 1)))]
        }
    } # end of deleting comments
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    ## Improved detection of the "Translate" command (2025-03-16) (the
    ## command must be preceded by at least one space OR a newline AND
    ## followed by the same)
    i2 <- grep("([[:space:]]+|^)TRANSLATE([[:space:]]+|$)",
               X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1) TRUE else FALSE
    if (length(i2) > 1)
        warning("problem translating labels in NEXUS file: is a taxon labelled \"Translate\"?")
    if (translation) {
        ## assume there's no semicolon in the labels:
        end <- semico[semico > i2][1]
        ## the old command:
        ## x <- X[(i2 + 1):end] # assumes there's a 'newline' after "TRANSLATE"
        ## the new command (doesn't assume TRANSLATE is on its own line):
        TRANS <- .decodeTRANSLATE(X[i2:end])
        n <- dim(TRANS)[1]
    }
    start <- if (translation) semico[semico > i2][1] + 1 else i1 + 1
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    rm(X)

    ## check whether there are empty lines from the above manips:
    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico)
    ## are some trees on several lines?
    ## -- this 'packs' all characters ending with a ";" in a single string
    if (Ntree == 1 && length(tree) > 1) {
        STRING <- paste(tree, collapse = "")
    } else {
        if (any(diff(semico) != 1)) {
            STRING <- character(Ntree)
            s <- c(1, semico[-Ntree] + 1)
            j <- mapply(":", s, semico)
            if (is.list(j)) {
                for (i in 1:Ntree)
                    STRING[i] <- paste(tree[j[[i]]], collapse = "")
            } else {
                for (i in 1:Ntree)
                    STRING[i] <- paste(tree[j[, i]], collapse = "")
            }
        } else STRING <- tree
    }
    rm(tree)
    ## exclude the possible command lines ending with ";":
    STRING <- STRING[grep("^[[:blank:]]*tree.*= *", STRING, ignore.case = TRUE)]
    Ntree <- length(STRING) # update Ntree
    if (is.null(tree.names)) {
        ## get the tree names:
        nms.trees <- sub(" *= *.*", "", STRING) # only the first occurence of "="
        nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", nms.trees, ignore.case = TRUE) # fix by Graham Gower (2014-10-20)
    }
    STRING <- sub("^.*= *", "", STRING) # delete title and 'TREE' command with 'sub'
    STRING <- gsub(" ", "", STRING) # delete all white spaces
    colon <- grep(":", STRING)

    FUN <- NULL
    if (length(colon) == 0) FUN <- .cladoBuild
    if (length(colon) == Ntree) FUN <- .treeBuild

    if (is.null(FUN)) {
        trees <- vector("list", Ntree)
        trees[colon] <- lapply(STRING[colon], .treeBuild)
        nocolon <- (1:Ntree)[-colon]
        trees[nocolon] <- lapply(STRING[nocolon], .cladoBuild)
    } else {
        trees <- lapply(STRING, FUN)
    }

    REF <- NULL
    if (translation) {
        for (i in 1:Ntree) {
            labs <- trees[[i]]$tip.label
            trees[[i]]$tip.label <- .translateTOKENS(labs, TRANS)
            labs <- trees[[i]]$node.label
            if (!is.null(labs))
                trees[[i]]$node.label <- .translateTOKENS(labs, TRANS)
        }
        if (identical(as.integer(TRANS[, 1]), 1:n))
            REF <- TRANS[, 2] # see GH's issue #141
    }

    if (Ntree == 1 && !force.multi) {
        if (!is.null(REF)) {
            tmp <- try(.compressTipLabel(trees, REF), silent = TRUE)
            if (inherits(tmp, "multiPhylo")) trees <- tmp
        }
        return(trees[[1]])
    }

    if (is.null(tree.names)) {
        if (!all(nms.trees == ""))
            names(trees) <- nms.trees
    } else {
        names(trees) <- tree.names
    }
    class(trees) <- "multiPhylo"
    tmp <- try(.compressTipLabel(trees, REF), silent = TRUE)
    if (inherits(tmp, "multiPhylo")) return(tmp)
    trees
}
