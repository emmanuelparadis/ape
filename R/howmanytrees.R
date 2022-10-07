## howmanytrees.R (2022-10-07)

##   Calculate Numbers of Phylogenetic Trees

## Copyright 2004-2022 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

LargeNumber <- function(a, b)
{
    c <- b * log10(a)
    n <- floor(c)
    x <- 10^(c - n)
    structure(c(x = x, n = n), class = "LargeNumber")
}

print.LargeNumber <- function(x, latex = FALSE, digits = 1, ...)
{
    if (latex) {
        cat("$", a, "^{", b, "} \\approx ", round(x, digits),
            " \\times 10^{", n, "}$\n", sep = "")
    } else {
        cat("approximately ", x["x"], " * 10^", x["n"], "\n", sep = "")
    }
}

howmanytrees <- function(n, rooted = TRUE, binary = TRUE,
                         labeled = TRUE, detail = FALSE)
{
    if (!labeled && !(rooted & binary))
      stop("can compute number of unlabeled trees only for rooted binary cases.")
    if (n < 3) N <- 1 else {
        if (labeled) {
            if (!rooted) n <- n - 1
            if (binary) {
                if (n < 152) {
                    N <- prod(seq(1, 2*n - 3, by = 2)) # double factorial
                } else {
                    ldfac <- lfactorial(2 * n - 3) - (n - 2) * log(2) - lfactorial(n - 2)
                    N <- LargeNumber(exp(1), ldfac)
                }
            }
            else {
                N <- matrix(0, n, n - 1)
                N[1:n, 1] <- 1
                for (i in 3:n)
                  for (j in 2:(i - 1))
                    N[i, j] <- (i + j - 2)*N[i - 1, j - 1] + j*N[i - 1, j]
                if (detail) {
                    rownames(N) <- 1:n
                    colnames(N) <- 1:(n - 1)
                } else N <- sum(N[n, ])
            }
        } else {
            N <- numeric(n)
            N[1] <- 1
            for (i in 2:n) {
                if (i %% 2) {
                    im1 <- i - 1L
                    x <- N[1:(im1 / 2)]
                    y <- N[im1:((i + 1) / 2)]
                } else {
                    ion2 <- i / 2
                    x <- N[1:ion2]
                    y <- N[(i - 1):ion2]
                    ny <- length(y)
                    y[ny] <- (y[ny] + 1) / 2
                }
                N[i] <- sum(x * y)
            }
            if (detail) names(N) <- 1:n else N <- N[n]
        }
    }
    N
}
