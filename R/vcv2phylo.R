## vcv2phylo.R (2014-11-27)

##   Variance-Covariance Matrix to Tree

## Copyright 2014 Simon Blomberg

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

vcv2phylo <- function (mat, tolerance = 1e-7)
{
  #########################################################
  ## Program to reconstruct a phylogenetic tree          ##
  ## from a phylogenetic variance-covariance matrix.     ##
  ## Input: mat (is tested for positive-definiteness)    ##
  ## Output: phylo (in phylo format as in package "ape") ##
  ## If numerical issues occur, adjust the tolerance.    ##
  ## Author: S. P. Blomberg                              ##
  ## Date: 12th November 2010                            ##
  #########################################################

  make.node <- function (left, right, value, lbrlen, rbrlen) {
    # function to make a node, using lists
    the.node <- list(left=left, right=right, value=value, lbrlen=lbrlen,
                     rbrlen=rbrlen)
    class(the.node) <- c("node", "list")
    return(the.node)
  }

  divide.matrix <- function (mat) {
    # function to decompose a block-diagonal matrix into
    # upper and lower blocks

    dims <- dim(mat)[1]
    end.of.block <- which(mat[1,] < tolerance)[1]-1
    if (is.na(end.of.block)) stop("Matrix is not block-diagonal")
    matlist <- list(upper=mat[1:end.of.block, 1:end.of.block],
                    lower=mat[(end.of.block+1):dims,(end.of.block+1):dims])
    if (length(matlist$upper)==1) names(matlist$upper) <- rownames(mat)[1]
    if (length(matlist$lower)==1) names(matlist$lower) <- rownames(mat)[dims]
    return(matlist)
  }

  make.tree.rec <- function (mat) {
    # Recursive function to create a tree made of nodes
    # from a phylogenetic matrix

    matlist <- divide.matrix(mat)
    if (is.vector(matlist$upper) && is.vector(matlist$lower)) {
      left <- as.numeric(names(matlist$upper))
      right <- as.numeric(names(matlist$lower))
      value <- i
      lbrlen <- matlist$upper
      rbrlen <- matlist$lower
    }
    if (is.vector(matlist$upper) && is.matrix(matlist$lower)) {
      min.lower <- min(matlist$lower)
      left <- as.numeric(names(matlist$upper))
      value <- i
      i <<- i+1
      right <- Recall(matlist$lower-min.lower)
      lbrlen <- matlist$upper
      rbrlen <- min.lower
    }
    if (is.matrix(matlist$upper) && is.vector(matlist$lower)) {
      min.upper <- min(matlist$upper)
      value <- i
      i <<- i+1
      left <- Recall(matlist$upper-min.upper)
      right <- as.numeric(names(matlist$lower))
      lbrlen <- min.upper
      rbrlen <- matlist$lower
    }
    if (is.matrix(matlist$upper) && is.matrix(matlist$lower)) {
      min.upper <- min(matlist$upper)
      min.lower <- min(matlist$lower)
      value <- i
      i <<- i+1
      left <- Recall(matlist$upper-min.upper)
      i <<- i+1
      right <- Recall(matlist$lower-min.lower)
      lbrlen <- min.upper
      rbrlen <- min.lower
    }
    return(make.node(left, right, value, lbrlen, rbrlen))
  }

  make.phylo.rec <- function (the.list) {
    # Recursive function to construct the edge matrix and collect the
    # branch length information from the tree
    brlens <<- c(brlens, the.list$lbrlen, the.list$rbrlen)
    if (is.numeric(the.list$left) && is.numeric(the.list$right)) {
      the.matrix <<- rbind(the.matrix, c(the.list$value, the.list$left),
                           c(the.list$value, the.list$right))
    }
    if (is.numeric(the.list$left) && inherits(the.list$right, "node")) {
      the.matrix <<- rbind(the.matrix, c(the.list$value, the.list$left),
                           c(the.list$value, the.list$right$value))
      Recall(the.list$right)
    }
    if (inherits(the.list$left, "node") && is.numeric(the.list$right)) {
        the.matrix <<- rbind(the.matrix, c(the.list$value, the.list$left$value),
                             c(the.list$value, the.list$right))
        Recall(the.list$left)
      }
    if (inherits(the.list$left, "node") && inherits(the.list$right, "node")) {
      the.matrix <<- rbind(the.matrix, c(the.list$value, the.list$left$value),
                           c(the.list$value, the.list$right$value))
      Recall(the.list$left)
      Recall(the.list$right)
    }
  }

  # main body
  #require(matrixcalc)
  #if (!is.positive.definite(mat)) stop("Matrix is not positive-definite")
  if (!isSymmetric(mat)) stop("Matrix is not symmetric")
  if (any(eigen(mat, only.values = TRUE)$values < -tolerance))
    stop("Matrix is not positive-definite")
  sp.names <- rownames(mat)
  dims <- dim(mat)[1]
  rownames(mat) <- colnames(mat) <- 1:dims
  i <- dims+1
  the.list <- make.tree.rec(mat) # side effect: calculate i
  the.matrix <- matrix(NA, 0, ncol=2) # initialise the edge matrix
  brlens <- vector(mode="numeric", length=0) #initialise branch length vector
  make.phylo.rec(the.list) # side effects: calculate the.matrix and brlens
  names(brlens) <- NULL
  phylo <- list(edge=the.matrix, tip.label=sp.names, edge.length=brlens, Nnode=i-dims)
  storage.mode(phylo$edge) <- "integer"
  storage.mode(phylo$Nnode) <- "integer"
  class(phylo) <- "phylo"
  return(phylo)
}

