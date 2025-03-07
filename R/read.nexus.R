## read.nexus.R (2024-07-22)

##   Read Tree File in Nexus Format

## Copyright 2003-2024 Emmanuel Paradis and 2010-2017 Klaus Schliep

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

.treeBuildWithTokens <- function(x)
{
    phy <- .Call(treeBuildWithTokens, x)
    dim(phy[[1]]) <- c(length(phy[[1]])/2, 2)
    nms <- c("edge", "edge.length", "Nnode", "node.label", "root.edge")
    if (length(phy) == 4) nms <- nms[-5]
    names(phy) <- nms
    if (all(phy$node.label == "")) phy$node.label <- NULL
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

## for read.nexus clado with TRANSLATION
.cladoBuildWithTokens <- function(x)
{
    phy <- .Call(cladoBuildWithTokens, x)
    dim(phy[[1]]) <- c(length(phy[[1]])/2, 2)
    nms <- c("edge", "Nnode", "node.label", "root.edge")
    if (length(phy) == 3) nms <- nms[-4]
    names(phy) <- nms
    if (all(phy$node.label == "")) phy$node.label <- NULL
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

.treeBuild <- function(x)
{
    if (!length(grep(",", x))) {
        ## only one tip but can be several nodes (GH's issues #104 and #124)
        Nnode <- length(gregexpr("\\)", x)[[1]])
        edge <- if (Nnode == 1L) 2:1 else c(2L, rep(3:(Nnode + 1L), each = 2), 1L)
        edge <- matrix(edge, Nnode, 2L, TRUE)
        phy <- list(edge = edge, Nnode = Nnode)
        labs <- unlist(strsplit(x, "[\\(\\):;]"))
        nt <- Nnode + 1L
        phy$tip.label <- labs[nt]
        labs <- labs[-(1:nt)]
        s <- c(TRUE, FALSE)
        tmp <- as.numeric(labs[s])
        if (length(tmp) == Nnode) {
            phy$edge.length <- tmp
        } else { # length(tmp) == Nnode + 1L (not checked)
            phy$edge.length <- tmp[-length(tmp)]
            phy$root.edge <- tmp[length(tmp)]
        }
        phy$node.label <- labs[!s]
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

read.nexus <- function(file, tree.names = NULL, force.multi = FALSE)
{
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)
    ## remove all comments
    ## (this might not work if there are square brackets within the comments)
    LEFT <- grep("\\[", X)
    RIGHT <- grep("\\]", X)
    if (length(LEFT)) { # in case there are no comments at all
        w <- LEFT == RIGHT
        if (any(w)) { # in case all comments use at least 2 lines
            s <- LEFT[w]
            X[s] <- gsub("\\[[^]]*\\]", "", X[s])
            ## The above regexp was quite tough to find: it makes
            ## possible to delete series of comments on the same line:
            ##       ...[...]xxx[...]...
            ## without deleting the "xxx". This regexp is in three parts:
            ##       \\[      [^]]*       \\]
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
    }
    endblock <- grep("END;|ENDBLOCK;", X, ignore.case = TRUE)
    semico <- grep(";", X)
    i1 <- grep("BEGIN TREES;", X, ignore.case = TRUE)
    i2 <- grep("TRANSLATE", X, ignore.case = TRUE)
    translation <- if (length(i2) == 1 && i2 > i1) TRUE else FALSE
    if (translation) {
        end <- semico[semico > i2][1]
        x <- X[(i2 + 1):end] # assumes there's a 'new line' after "TRANSLATE"
        #x <- unlist(strsplit(x, "[,; \t]"))
        ################################################
        # when the label of translation contains space #
        # 1 "tip 1 a",                                 #
        # 2 "tip 2"                                    #
        ################################################
        # remove the space and tab before the string
        x <- gsub("^\\s+", "", x)
        # remove the , ; symbol
        x <- gsub("[,;]", "", x)
        # split with the first space
        x <- unlist(regmatches(x, regexpr("\\s+", x), invert=TRUE))
        ###############################################
        x <- x[nzchar(x)]
        TRANS <- matrix(x, ncol = 2, byrow = TRUE)
        TRANS[, 2] <- gsub("['\"]", "", TRANS[, 2])
        n <- dim(TRANS)[1]
    }
    trans_vec <- setNames(TRANS[,2], TRANS[,1])
    start <-
        if (translation) semico[semico > i2][1] + 1
        else i1 + 1 # semico[semico > i1][1] ## fix done on 2014-08-25
    end <- endblock[endblock > i1][1] - 1
    tree <- X[start:end]
    rm(X)

    ## check whether there are empty lines from the above manips:
    tree <- tree[tree != ""]
    semico <- grep(";", tree)
    Ntree <- length(semico) # provisional -- some ";" may actually mark end of commands
    ## are some trees on several lines?
    ## -- this actually 'packs' all characters ending with a ";" in a single string --
    if (Ntree == 1 && length(tree) > 1) STRING <- paste(tree, collapse = "") else {
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
    ## get the tree names:
    nms.trees <- sub(" *= *.*", "", STRING) # only the first occurence of "="
    nms.trees <- sub("^[[:blank:]]*tree[[:blank:]\\*]*", "", nms.trees, ignore.case = TRUE) # fix by Graham Gower (2014-10-20)
    STRING <- sub("^.*= *", "", STRING) # delete title and 'TREE' command with 'sub'
    STRING <- gsub(" ", "", STRING) # delete all white spaces
    colon <- grep(":", STRING)
    
    with_token <- FALSE
    
    if (!length(colon)) {
        trees <- lapply(STRING, .cladoBuild)
    } else if (length(colon) == Ntree) {
        if (translation) with_token <- all(lengths(gregexpr("\\,", STRING)) == (length(trans_vec)-1L) )
        trees <-
            if (with_token) lapply(STRING, .treeBuildWithTokens)
            else lapply(STRING, .treeBuild)
    } else {
        trees <- vector("list", Ntree)
        trees[colon] <- lapply(STRING[colon], .treeBuild)
        nocolon <- (1:Ntree)[!1:Ntree %in% colon]
        trees[nocolon] <- lapply(STRING[nocolon], .cladoBuild)
        if (translation) {
            for (i in 1:Ntree) {
                tr <- trees[[i]]
                for (j in 1:n) {
                    ind <- which(tr$tip.label[j] == TRANS[, 1])
                    tr$tip.label[j] <- TRANS[ind, 2]
                }
                if (!is.null(tr$node.label)) {
                    for (j in 1:length(tr$node.label)) {
                        ind <- which(tr$node.label[j] == TRANS[, 1])
                        tr$node.label[j] <- TRANS[ind, 2]
                    }
                }
                trees[[i]] <- tr
            }
            translation <- FALSE
        }
    }
    for (i in 1:Ntree) {
        tr <- trees[[i]]
        if (!translation) n <- length(tr$tip.label)
    }
    if (Ntree == 1 && !force.multi) {
        trees <- trees[[1]]
        if (translation) {
          trees$tip.label <- as.vector(trans_vec[trees$tip.label])
#            trees$tip.label <-
#                if (length(colon)) TRANS[, 2] else
#                TRANS[, 2][as.numeric(trees$tip.label)]
        }
    } else {
        if (!is.null(tree.names)) names(trees) <- tree.names
        if (translation) {
            if (with_token) # && length(colon) == Ntree) # .treeBuildWithTokens() was used
                attr(trees, "TipLabel") <- TRANS[, 2]
            else { # reassign the tip labels then compress
                for (i in 1:Ntree)
                    trees[[i]]$tip.label <- as.vector(trans_vec[trees[[i]]$tip.label]) 
#                    trees[[i]]$tip.label <-
#                        TRANS[, 2][as.numeric(trees[[i]]$tip.label)]
                class(trees) <- "multiPhylo"
                if(all(Ntip(trees)==length(trans_vec))) trees <- .compressTipLabel(trees)
            }
        }
        class(trees) <- "multiPhylo"
        if (!all(nms.trees == "")) names(trees) <- nms.trees
    }
    trees
}
