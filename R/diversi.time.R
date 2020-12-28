## diversi.time.R (2007-09-22)

##   Analysis of Diversification with Survival Models

## Copyright 2002-2007 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

diversi.time <- function(x, census = NULL, censoring.codes = c(1, 0),
                         Tc = NULL)
{
    n <- length(x)
    if (is.null(census)) {
        k <- n
        census <- rep(censoring.codes[1], n)
    }
    else k <- sum(census == censoring.codes[1])
    u <- n - k
    S <- sum(x)
    delta <- k / S
    var.delta <- delta^2 / k
    loglik.A <- k * log(delta) - delta * S
    tk <- x[census == censoring.codes[1]]
    tu <- x[census == censoring.codes[2]]
    fb <- function(b)
      1/b - sum(x^b * log(x))/sum(x^b) + sum(log(tk))/k
    beta <- uniroot(fb, interval = c(1e-7, 10))$root
    Sp <- sum(x^beta)
    alpha <- (k / Sp)^(1/beta)
    var.alpha <- 1/ ((k * beta / alpha^2) + beta * (beta - 1) * alpha^(beta - 2) * Sp)
    ax <- alpha * x
    var.beta <- 1 / (k / beta^2 + sum(ax^beta * log(ax)))
    loglik.B <- k*(log(alpha) + log(beta)) +
      (beta - 1)*(k*log(alpha) + sum(log(tk)))- Sp * alpha^beta
    if (is.null(Tc)) Tc <- median(x)
    tk1 <- tk[tk < Tc]
    tk2 <- tk[tk >= Tc]
    tu1 <- tu[tu < Tc]
    tu2 <- tu[tu >= Tc]
    k1 <- length(tk1)
    k2 <- k - k1
    u1 <- length(tu1)
    u2 <- u - u1
    tmp <- (k2 + u2) * Tc
    delta1 <- k1 / (sum(tk1) + sum(tu1) + tmp)
    delta2 <- k2 / (sum(tk2) + sum(tu2) - tmp)
    var.delta1 <- delta1^2 / k1
    var.delta2 <- delta2^2 / k2
    tmp <- Tc * (delta2 - delta1)
    loglik.C <- k1 * log(delta1) - delta1 * sum(tk1) + k2 * log(delta2) +
                  k2 * tmp - delta2 * sum(tk2) - delta1 * sum(tu1) +
                    u2 * tmp - delta2 * sum(tu2)
    cat("\nAnalysis of Diversification with Survival Models\n\n")
    cat("Data:", deparse(substitute(x)), "\n")
    cat("Number of branching times:", n, "\n")
    cat("         accurately known:", k, "\n")
    cat("                 censored:", u, "\n\n")
    cat("Model A: constant diversification\n")
    cat("    log-likelihood =", round(loglik.A, 3),
        "   AIC =", round(-2 * loglik.A + 2, 3), "\n")
    cat("    delta =", round(delta, 6), "   StdErr =", round(sqrt(var.delta), 6), "\n\n")
    cat("Model B: diversification follows a Weibull law\n")
    cat("    log-likelihood =", round(loglik.B, 3),
        "   AIC =", round(-2 * loglik.B + 4, 3), "\n")
    cat("    alpha =", round(alpha, 6), "   StdErr =", round(sqrt(var.alpha), 6), "\n")
    cat("     beta =", round(beta, 6), "   StdErr =", round(sqrt(var.beta), 6), "\n\n")
    cat("Model C: diversification changes with a breakpoint at time =", Tc, "\n")
    cat("    log-likelihood =", round(loglik.C, 3),
        "   AIC =", round(-2 * loglik.C + 4, 3), "\n")
    cat("    delta1 =", round(delta1, 6), "   StdErr =", round(sqrt(var.delta1), 6), "\n")
    cat("    delta2 =", round(delta2, 6), "   StdErr =", round(sqrt(var.delta2), 6), "\n\n")
    cat("Likelihood ratio tests:\n")
    c1 <- 2 * (loglik.B - loglik.A)
    p1 <- round(1 - pchisq(c1, 1), 4)
    c2 <- 2 * (loglik.C - loglik.A)
    p2 <- round(1 - pchisq(c2, 1), 4)
    cat("    Model A vs. Model B: chi^2 =", round(c1, 3), "   df = 1,    P =", p1, "\n")
    cat("    Model A vs. Model C: chi^2 =", round(c2, 3), "   df = 1,    P =", p2, "\n")
}
