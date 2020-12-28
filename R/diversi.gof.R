## diversi.gof.R (2006-10-16)

##   Tests of Constant Diversification Rates

## Copyright 2002-2006 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

diversi.gof <- function(x, null = "exponential", z = NULL)
{
    n <- length(x)
    if (null == "exponential") {
        delta <- n/sum(x)
        z <- 1 - exp(-delta * sort(x))
    }
    else {
        nmsz <- deparse(substitute(z))
        z <- sort(z) # utile ???
    }
    i <- 1:n
    W2 <- sum((z - (2*i - 1)/(2*n))^2) + 1/12*n
    A2 <- -sum((2*i - 1)*(log(z) + log(1 - rev(z))))/n - n
    if (null == "exponential") {
        W2 <- W2*(1 - 0.16/n)
        A2 <- A2*(1 + 0.6/n)
    }
    else W2 <- (W2 - 0.4/n + 0.6/n^2)/(1 + 1/n)
    cat("\nTests of Constant Diversification Rates\n\n")
    cat("Data:", deparse(substitute(x)), "\n")
    cat("Number of branching times:", n, "\n")
    cat("Null model: ")
    if (null == "exponential") cat("exponential\n\n")
    else cat(nmsz, "(user-specified)\n\n")
    cat("Cramer-von Mises test: W2 =", round(W2, 3))
    if (null == "exponential") {
        if (W2 < 0.177) cat("   P > 0.1\n")
        if (W2 >= 0.177 && W2 < 0.224) cat("   0.05 < P < 0.1\n")
        if (W2 >= 0.224 && W2 < 0.273) cat("   0.025 < P < 0.05\n")
        if (W2 >= 0.273 && W2 < 0.337) cat("   0.01 < P < 0.025\n")
        if (W2 > 0.337) cat("   P < 0.01\n")
    }
    else {
        if (W2 < 0.347) cat("   P > 0.1\n")
        if (W2 >= 0.347 && W2 < 0.461) cat("   0.05 < P < 0.1\n")
        if (W2 >= 0.461 && W2 < 0.581) cat("   0.025 < P < 0.05\n")
        if (W2 >= 0.581 && W2 < 0.743) cat("   0.01 < P < 0.025\n")
        if (W2 > 0.743) cat("   P < 0.01\n")
    }
    cat("Anderson-Darling test: A2 =", round(A2, 3))
    if (null == "exponential") {
        if (A2 < 1.078) cat("   P > 0.1\n")
        if (A2 >= 1.078 && A2 < 1.341) cat("   0.05 < P < 0.1\n")
        if (A2 >= 1.341 && A2 < 1.606) cat("   0.025 < P < 0.05\n")
        if (A2 >= 1.606 && A2 < 1.957) cat("   0.01 < P < 0.025\n")
        if (A2 > 1.957) cat("   P < 0.01\n")
    }
    else {
        if (A2 < 1.933) cat("   P > 0.1\n")
        if (A2 >= 1.933 && A2 < 2.492) cat("   0.05 < P < 0.1\n")
        if (A2 >= 2.492 && A2 < 3.070) cat("   0.025 < P < 0.05\n")
        if (A2 >= 3.070 && A2 < 3.857) cat("   0.01 < P < 0.025\n")
        if (A2 > 3.857) cat("   P < 0.01\n")
    }
}
