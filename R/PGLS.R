## PGLS.R (2020-04-13)

##   Phylogenetic Generalized Least Squares

## Copyright 2004-2020 Julien Dutheil, and 2006-2017 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

corBrownian <- function(value = 1, phy, form = ~1)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- TRUE
    attr(value, "tree") <- phy
    class(value) <- c("corBrownian", "corPhyl", "corStruct")
    value
}

corMartins <- function(value, phy, form = ~1, fixed = FALSE)
{
    if (length(value) > 1)
        stop("only one parameter is allowed")
    if (value < 0) stop("the parameter alpha must be positive")
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    attr(value, "tree") <- phy
    class(value) <- c("corMartins", "corPhyl", "corStruct")
    value
}

corGrafen <- function(value, phy, form = ~1, fixed = FALSE)
{
    if (length(value) > 1)
        stop("only one parameter is allowed")
    if (value < 0) stop("parameter rho must be positive")
    value <- log(value) # Optimization under constraint, use exponential transform.
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    attr(value, "tree") <- phy
    class(value) <- c("corGrafen", "corPhyl", "corStruct")

    value
}

Initialize.corPhyl <- function(object, data, ...)
{
    ## The same as in Initialize corStruct:
    form <- formula(object)
    if (getCovariateFormula(object) == ~1) {
      warning("No covariate specified, species will be taken as ordered in the data frame. To avoid this message, specify a covariate containing the species names with the 'form' argument.")
    }
    ## Obtaining the group information, if any
    if (!is.null(getGroupsFormula(form))) {
        attr(object, "groups") <- getGroups(object, form, data = data)
        attr(object, "Dim") <- Dim(object, attr(object, "groups"))
    } else { # no groups
        attr(object, "Dim") <- Dim(object, as.factor(rep(1, nrow(data))))
    }
    ## Obtaining the covariate(s)
    attr(object, "covariate") <- getCovariate(object, data = data)

    object
}

corMatrix.corBrownian <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (!("corBrownian" %in% class(object)))
        stop('object is not of class "corBrownian"')
    if (data.class(covariate) == "list") {
        as.list(lapply(covariate, function(el) corMatrix(object, covariate = el)))
    } else {
        tree <- attr(object, "tree")
        mat <- vcv.phylo(tree, corr = corr)
        mat[covariate, covariate]
    }
}

corMatrix.corMartins <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (!("corMartins" %in% class(object)))
        stop('object is not of class "corMartins"')
    if (data.class(covariate) == "list") {
        as.list(lapply(covariate, function(el) corMatrix(object, covariate = el)))
    } else {
        tree <- attr(object, "tree")
        dist <- cophenetic.phylo(tree)
        mat <- exp(-object[1] * dist)
        if (corr) mat <- cov2cor(mat)
        mat[covariate, covariate]
    }
}

corMatrix.corGrafen <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (!("corGrafen" %in% class(object)))
        stop('object is not of class "corGrafen"')
    if (data.class(covariate) == "list") {
        as.list(lapply(covariate, function(el) corMatrix(object, covariate = el)))
    } else {
        tree <- compute.brlen(attr(object, "tree"),
                              method = "Grafen", power = exp(object[1]))
        mat <- vcv.phylo(tree, corr = corr)
        mat[covariate, covariate]
    }
}

coef.corBrownian <- function(object, unconstrained = TRUE, ...)
{
    if (!("corBrownian" %in% class(object)))
        stop('object is not of class "corBrownian"')
    numeric(0)
}

coef.corMartins <- function(object, unconstrained = TRUE, ...)
{
    if (!("corMartins" %in% class(object)))
        stop('object is not of class "corMartins"')
    if (unconstrained) {
        if (attr(object, "fixed")) {
            return(numeric(0))
        } else {
            return(as.vector(object))
        }
    }
    aux <- as.vector(object)
    names(aux) <- "alpha"
    aux
}

coef.corGrafen <- function(object, unconstrained = TRUE, ...)
{
    if (!("corGrafen" %in% class(object)))
        stop('object is not of class "corGrafen"')
    if (unconstrained) {
        if (attr(object, "fixed")) {
            return(numeric(0))
        } else {
            return(as.vector(object))
        }
    }
    aux <- exp(as.vector(object))
    names(aux) <- "rho"
    aux
}

## changed by EP (2006-10-12):

compute.brlen <- function(phy, method = "Grafen", power = 1, ...)
{
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    Ntip <- length(phy$tip.label)
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (is.numeric(method)) {
        phy$edge.length <- rep(method, length.out = Nedge)
        return(phy)
    }
    if (is.function(method)) {
        phy$edge.length <- method(Nedge, ...)
        return(phy)
    }
    if (is.character(method)) { # == "Grafen"
        tr <- reorder(phy, "postorder")
        xx <- .C(node_depth, as.integer(Ntip),
                 as.integer(tr$edge[, 1]), as.integer(tr$edge[, 2]),
                 as.integer(Nedge), double(Ntip + Nnode), 1L)[[5]] - 1
        m <- Ntip - 1
        phy$edge.length <-
          (xx[phy$edge[, 1]]/m)^power - (xx[phy$edge[, 2]]/m)^power
        return(phy)
    }
}

## by EP:

corPagel <- function(value, phy, form = ~1, fixed = FALSE)
{
    if (value < 0 || value > 1)
        stop("the value of lambda must be between 0 and 1.")
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    attr(value, "tree") <- phy
    class(value) <- c("corPagel", "corPhyl", "corStruct")
    value
}

corMatrix.corPagel <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (!("corPagel" %in% class(object)))
        stop('object is not of class "corPagel"')
    if (data.class(covariate) == "list") {
        as.list(lapply(covariate, function(el) corMatrix(object, covariate = el)))
    } else {
        mat <- vcv.phylo(attr(object, "tree"), corr = corr)
        mat <- mat[covariate, covariate]
        tmp <- diag(mat)
        mat <- object[1]*mat
        diag(mat) <- tmp
        mat
    }
}

coef.corPagel <- function(object, unconstrained = TRUE, ...)
{
    if (unconstrained) {
        if (attr(object, "fixed")) return(numeric(0))
        else return(object[1])
    }
    aux <- object[1]
    names(aux) <- "lambda"
    aux
}

corBlomberg <- function(value, phy, form = ~1, fixed = FALSE)
{
    if (value <= 0)
        stop("the value of g must be greater than 0.")
    if (!inherits(phy, "phylo"))
        stop('object "phy" is not of class "phylo"')
    attr(value, "formula") <- form
    attr(value, "fixed") <- fixed
    attr(value, "tree") <- phy
    class(value) <- c("corBlomberg", "corPhyl", "corStruct")
    value
}

corMatrix.corBlomberg <-
    function(object, covariate = getCovariate(object), corr = TRUE, ...)
{
    if (object[1] <= 0)
        stop("the optimization has reached a value <= 0 for parameter 'g':
probably need to set 'fixed = TRUE' in corBlomberg().")
    if (data.class(covariate) == "list") {
        as.list(lapply(covariate, function(el) corMatrix(object, covariate = el)))
    } else {
        phy <- attr(object, "tree")
        d <- (dist.nodes(phy)[length(phy$tip.label) + 1, ])^(1/object[1])
        phy$edge.length <- d[phy$edge[, 2]] - d[phy$edge[, 1]]
        mat <- vcv.phylo(phy, corr = corr)
        mat[covariate, covariate]
    }
}

coef.corBlomberg <- function(object, unconstrained = TRUE, ...)
{
    if (unconstrained) {
        if (attr(object, "fixed")) return(numeric(0))
        else return(object[1])
    }
    aux <- object[1]
    names(aux) <- "g"
    aux
}
