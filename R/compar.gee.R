## compar.gee.R (2015-05-01)

##   Comparative Analysis with GEEs

## Copyright 2002-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

compar.gee <-
    function(formula, data = NULL, family = gaussian, phy,
             corStruct, scale.fix = FALSE, scale.value = 1)
{
    if (requireNamespace("gee", quietly = TRUE)) gee <- gee::gee
    else stop("package 'gee' not available")

    if (!missing(corStruct)) {
        if (!missing(phy))
            warning("the phylogeny was ignored because you gave a 'corStruct' object")
        R <- vcv(corStruct, corr = TRUE)
    } else {
        R <- vcv(phy, corr = TRUE)
    }

    if (is.null(data)) data <- parent.frame()
    else {
        nmsR <- rownames(R)
        if (!identical(rownames(data), nmsR)) {
            if (!any(is.na(match(rownames(data), nmsR))))
                data <- data[nmsR, ]
            else {
                msg <- if (missing(corStruct))
                    "the tip labels of the tree" else "those of the correlation structure"
                msg <- paste("the rownames of the data.frame and", msg,
                             "do not match: the former were ignored in the analysis")
                warning(msg)
            }
        }
    }

    effect.assign <- attr(model.matrix(formula, data = data), "assign")

    for (i in all.vars(formula)) {
        if (any(is.na(eval(parse(text = i), envir = data))))
          stop("the present method cannot be used with missing data: you may consider removing the species with missing data from your tree with the function 'drop.tip'.")
    }

    id <- rep(1, dim(R)[1])
    geemod <- do.call("gee", list(formula, id, data = data, family = family, R = R,
                                  corstr = "fixed", scale.fix = scale.fix,
                                  scale.value = scale.value))
    W <- geemod$naive.variance
    fname <-
        if (is.function(family)) deparse(substitute(family)) else if (is.list(family)) family$family else family
    if (fname == "binomial")
        W <- summary(glm(formula, family = quasibinomial, data = data))$cov.scaled
    N <- geemod$nobs
    ## <FIXME>
    ## maybe need to refine below in case of non-Brownian corStruct
    if (!missing(corStruct)) phy <- attr(corStruct, "tree")
    dfP <- sum(phy$edge.length)*N / sum(diag(vcv(phy))) # need the variances
    ## </FIXME>

    ## compute QIC:
    Y <- geemod$y
    MU <- geemod$fitted.values
    Qlik <- switch(fname,
                   "gaussian" = -sum((Y - MU)^2)/2,
                   "binomial" = sum(Y*log(MU/(1 - MU)) + log(1 - MU)),
                   "poisson" = sum(Y*log(MU) - MU),
                   "Gamma" = sum(Y/MU + log(MU)),
                   "inverse.gaussian" = sum(-Y/(2*MU^2) + 1/MU))
    Ai <- do.call("gee", list(formula, id, data = data, family = family,
                              corstr = "independence", scale.fix = scale.fix,
                              scale.value = scale.value))$naive.variance
    QIC <- -2*Qlik + 2*sum(diag(solve(Ai) %*% W))

    obj <- list(call = match.call(),
                effect.assign = effect.assign,
                nobs = N,
                QIC = QIC,
                coefficients = geemod$coefficients,
                residuals = geemod$residuals,
                fitted.values = MU,
                family = geemod$family$family,
                link = geemod$family$link,
                scale = geemod$scale,
                W = W,
                dfP = dfP)
    class(obj) <- "compar.gee"
    obj
}

print.compar.gee <- function(x, ...)
{
    nas <- is.na(x$coef)
    coef <- x$coef[!nas]
    cnames <- names(coef)
    coef <- matrix(rep(coef, 4), ncol = 4)
    dimnames(coef) <- list(cnames,
                           c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
    df <- x$dfP - dim(coef)[1]
    coef[, 2] <- sqrt(diag(x$W))
    coef[, 3] <- coef[, 1]/coef[, 2]
    if (df < 0) {
        warning("not enough degrees of freedom to compute P-values.")
        coef[, 4] <- NA
    } else coef[, 4] <- 2 * (1 -  pt(abs(coef[, 3]), df))
    residu <- quantile(as.vector(x$residuals))
    names(residu) <- c("Min", "1Q", "Median", "3Q", "Max")
    cat("Call: ")
    print(x$call)
    cat("Number of observations: ", x$nobs, "\n")
    cat("Model:\n")
    cat("                      Link:", x$link, "\n")
    cat(" Variance to Mean Relation:", x$family, "\n")
    cat("\nQIC:", x$QIC, "\n")
    cat("\nSummary of Residuals:\n")
    print(residu)
    if (any(nas))
        cat("\n\nCoefficients: (", sum(nas), " not defined because of singularities)\n",
            sep = "")
    else cat("\n\nCoefficients:\n")
    print(coef)
    cat("\nEstimated Scale Parameter: ", x$scale)
    cat("\n\"Phylogenetic\" df (dfP): ", x$dfP, "\n")
}

drop1.compar.gee <- function(object, scope, quiet = FALSE, ...)
{
    fm <- formula(object$call)
    trm <- terms(fm)
    z <- attr(trm, "term.labels")
    ind <- object$effect.assign
    n <- length(z)
    ans <- matrix(NA, n, 3)
    for (i in 1:n) {
        wh <- which(ind == i)
        ans[i, 1] <- length(wh)
        ans[i, 2] <- t(object$coefficients[wh]) %*%
          solve(object$W[wh, wh]) %*% object$coefficients[wh]
    }
    df <- object$dfP - length(object$coefficients)
    if (df < 0) warning("not enough degrees of freedom to compute P-values.")
    else ans[, 3] <- pf(ans[, 2], ans[, 1], df, lower.tail = FALSE)
    colnames(ans) <- c("df", "F", "Pr(>F)")
    rownames(ans) <- z
    if (any(attr(trm, "order") > 1) && !quiet)
      warning("there is at least one interaction term in your model:
you should be careful when interpreting the significance of the main effects.")
    class(ans) <- "anova"
    attr(ans, "heading") <- paste("Single term deletions\n\n  Model:",
                                  as.character(as.expression(fm)), "\n")
    ans
}

predict.compar.gee <-
    function(object, newdata = NULL, type = c("link", "response"), ...)
{
    type <- match.arg(type)
    pred <- if (is.null(newdata)) object$fitted.values else {
        frm <- formula(object$call$formula)[-2]
        X <-  model.matrix(frm, data = newdata)
        beta <- object$coefficients
        X[, names(beta), drop = FALSE] %*% beta
    }
    if (type == "link") return(pred)
    f <- match.fun(object$family)
    f(link = object$link)$linkinv(pred)
}
