## dbd.R (2015-02-06)

##   Probability Density Under Birth--Death Models

## Copyright 2012-2015 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

dyule <- function(x, lambda = 0.1, t = 1, log = FALSE)
{
    tmp <- exp(-lambda * t)
    res <- if (log) log(tmp) + (x - 1) * log(1 - tmp) else tmp * (1 - tmp)^(x - 1)
    out.of.range <- x <= 0
    if (any(out.of.range))
        res[out.of.range] <- if (log) -Inf else 0
    res
}

dbd <- function(x, lambda, mu, t, conditional = FALSE, log = FALSE)
{
    if (length(lambda) > 1) {
        lambda <- lambda[1]
        warning("only the first value of 'lambda' was considered")
    }
    if (length(mu) > 1) {
        mu <- mu[1]
        warning("only the first value of 'mu' was considered")
    }

    if (mu == 0) return(dyule(x, lambda, t, log))

    ## for the unconditional case, we have to consider x=0 separately:
    if (!conditional) {
        zero <- x == 0
        out.of.range <- x < 0
    } else {
        out.of.range <- x <= 0
    }

    res <- numeric(length(x))

    ## the situation were speciation and extinction probabilities are equal:
    if (lambda == mu) {
        tmp <- lambda * t
        eta <- tmp/(1 + tmp)
        if (conditional) {
            res[] <- if (log) log(1 - eta) + (x - 1) * log(eta) else (1 - eta) * eta^(x - 1)
        } else { # the unconditional case:
            if (length(zero)) {
                res[zero] <- eta
                res[!zero] <- (1 - eta)^2 * eta^(x[!zero] - 1)
            } else res[] <- (1 - eta)^2 * eta^(x - 1)
        }
    } else { # the general case with lambda != mu

        ## this expression is common to the conditional and unconditional cases:
        Ent <- exp((lambda - mu) * t)

        if (conditional) {
            if (log) {
                res[] <- log(lambda - mu) - log(lambda * Ent - mu) +
                    (x - 1) * (log(lambda) + log(Ent - 1) - log(lambda * Ent - mu))
            } else {
                eta <- lambda * (Ent - 1)/(lambda * Ent - mu)
                res[] <- (1 - eta) * eta^(x - 1)
            }
        } else { # finally, the unconditional case:
            eta <- lambda * (Ent - 1)/(lambda * Ent - mu)
            if (length(zero)) {
                res[zero] <- eta * mu / lambda
                res[!zero] <- (1 - mu * eta / lambda) * (1 - eta) * eta^(x[!zero] - 1)
            } else res[] <- (1 - mu * eta / lambda) * (1 - eta) * eta^(x - 1)
        }
    }

    if (any(out.of.range))
        res[out.of.range] <- if (log) -Inf else 0
    res
}

dbdTime <- function(x, birth, death, t, conditional = FALSE,
                    BIRTH = NULL, DEATH = NULL, fast = FALSE)
{
    if (length(t) > 1) {
        t <- t[1]
        warning("only the first value of 't' was considered")
    }

    if (conditional) {
        PrNt <- function(t, T, x) {
            tmp <- exp(-RHO(t, T))
            Wt <- tmp * (1 + INT(t))
            out <- (1/Wt)*(1 - 1/Wt)^(x - 1)
            zero <- x == 0
            if (length(zero)) out[zero] <- 0
            out
        }
    } else { # the unconditional case:
        PrNt <- function(t, T, x)
        {
            tmp <- exp(-RHO(t, T))
            Wt <- tmp * (1 + INT(t))
            out <- numeric(length(x))
            zero <- x == 0
            if (length(zero)) {
                out[zero] <- 1 - tmp/Wt
                out[!zero] <- (tmp/Wt^2)*(1 - 1/Wt)^(x[!zero] - 1)
            } else out[] <- (tmp/Wt^2)*(1 - 1/Wt)^(x - 1)
            out
        }
    }
    case <- .getCase(birth, death, BIRTH, DEATH)
    ff <- .getRHOetINT(birth, death, BIRTH, DEATH, case = case, fast = fast)
    RHO <- ff[[1]]
    INT <- ff[[2]]
    environment(RHO) <- environment(INT) <- environment()
    Tmax <- t
    PrNt(0, t, x)
}
