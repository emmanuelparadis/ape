## SDM.R (2012-04-02)

## Construction of Consensus Distance Matrix With SDM

## Copyright 2011-2012 Andrei-Alin Popescu

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

SDM <- function(...)
{
    st <- list(...) # first half contains matrices, second half s_p
    k <- length(st)/2
    ONEtoK <- seq_len(k)

    ## make sure we have only matrices:
    for (i in ONEtoK) st[[i]] <- as.matrix(st[[i]])

    ## store the rownames of each matrix in a list because they are often called:
    ROWNAMES <- lapply(st[ONEtoK], rownames)

    ## the number of rows of each matrix:
    NROWS <- lengths(ROWNAMES)
    tot <- sum(NROWS)

    labels <- unique(unlist(ROWNAMES))
    sp <- unlist(st[k + ONEtoK])

    astart <- numeric(tot) # start of aip, astart[p] is start of aip
    astart[1] <- k
    for (i in 2:k)
        astart[i] <- astart[i - 1] + NROWS[i - 1]

    ## apparently erased by the operation below so no need to initialize:
    ## a <- mat.or.vec(1, k + tot + k + length(labels))

    ## first k are alphas, subsequent ones aip
    ## each matrix p starting at astart[p], next are
    ## Lagrange multipliers, miu, niu, lambda in that order
    n <- length(labels)
    miustart <- k + tot
    niustart <- miustart + n
    lambstart <- niustart + k - 1

    X <- matrix(0, n, n, dimnames = list(labels, labels))
    V <- w <- X

    tmp <- 2 * k + tot + n
    col <- numeric(tmp) # free terms of system

    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            for (p in ONEtoK) {
                ## d <- st[[p]] # not needed anymore
                if (is.element(labels[i], ROWNAMES[[p]]) && is.element(labels[j], ROWNAMES[[p]])) {
                    w[i, j] <- w[j, i] <- w[i, j] + sp[p]
                }
            }
        }
    }

    ONEtoN <- seq_len(n)

    Q <- matrix(0, tmp, tmp)
    ## first decompose first sum in paper
    for (p in ONEtoK) {
        d_p <- st[[p]]
        for (l in ONEtoK) { # first compute coefficients of alphas
            d <- st[[l]]
            sum <- 0
            dijp <- -1
            if (l == p) { # calculate alpha_p
                for (i in ONEtoN) {
                    for (j in ONEtoN) { # check if {i,j}\subset L_l
                        if (i == j) next # make sure i != j
                        ## d <- st[[l]] # <- moved-up
                        pos <- match(labels[c(i, j)], ROWNAMES[[l]]) # <- returns NA if not in this matrix
                        if (all(!is.na(pos))) {
                            ipos <- pos[1L]
                            jpos <- pos[2L]
                            dij <- d[ipos, jpos]
                            sum <- sum + dij * dij - sp[l] * dij * dij / w[i,j]
                            tmp2 <- dij - sp[l] * dij / w[i,j]
                            Q[p, astart[l] + ipos] <- Q[p, astart[l] + ipos] + tmp2
                            Q[p, astart[l] + jpos] <- Q[p, astart[l] + jpos] + tmp2
                        }
                    }
                }
            } else {
                for (i in ONEtoN) {
                    for (j in ONEtoN) { # check if {i,j}\subset L_l
                        if (i == j) next
                        ## d <- st[[l]] # <- moved-up
                        pos <- match(labels[c(i, j)], ROWNAMES[[l]])
                        posp <- match(labels[c(i, j)], ROWNAMES[[p]])
                        if (all(!is.na(pos)) && all(!is.na(posp))) {
                            ipos <- pos[1L]
                            jpos <- pos[2L]
                            dij <- d[ipos, jpos]
                            dijp <- d_p[posp[1L], posp[2L]]
                            sum <- sum - sp[l] * dij * dijp / w[i, j]
                            tmp2 <- sp[l] * dijp / w[i, j]
                            Q[p,astart[l] + ipos] <- Q[p, astart[l] + ipos] - tmp2
                            Q[p,astart[l] + jpos] <- Q[p, astart[l] + jpos] - tmp2
                        }
                    }
                }
            }
            Q[p, l] <- sum
        }
        Q[p, lambstart + 1] <- 1
    }

    r <- k

    for (p in ONEtoK) {
        dp <- st[[p]]
        for (i in ONEtoN) {
            if (is.element(labels[i], ROWNAMES[[p]])) {
                r <- r + 1
                for (l in ONEtoK) {
                    d <- st[[l]]
                    if (l == p) {
                        ipos <- match(labels[i], ROWNAMES[[p]])
                        for (j in ONEtoN) {
                            if (i == j) next
                            jpos <- match(labels[j], ROWNAMES[[p]])
                            if (!is.na(jpos)) {
                                dij <- d[ipos, jpos]
                                Q[r, l] <- Q[r, l] + dij - sp[l] * dij / w[i, j]
                                tmp2 <- 1 - sp[l] / w[i, j]
                                Q[r, astart[l] + ipos] <- Q[r, astart[l] + ipos] + tmp2
                                Q[r, astart[l] + jpos] <- Q[r, astart[l] + jpos] + tmp2
                            }
                        }
                    } else {
                        for (j in ONEtoN) {
                            if (i == j) next
                            if (!is.element(labels[j], rownames(dp))) next
                            pos <- match(labels[c(i, j)], ROWNAMES[[l]])
                            if (all(!is.na(pos))) {
                                ipos <- pos[1L]
                                jpos <- pos[2L]
                                dij <- d[ipos, jpos]
                                Q[r, l] <- Q[r, l] - sp[l] * dij / w[i, j]
                                tmp2 <- sp[l]/w[i, j]
                                Q[r, astart[l] + ipos] <- Q[r, astart[l] + ipos] - tmp2
                                Q[r, astart[l] + jpos] <- Q[r, astart[l] + jpos] - tmp2
                            }
                        }
                    }
                }
                if (p < k) Q[r, ] <- Q[r, ] * sp[p]
                Q[r, miustart + i] <- 1
                if (p < k) Q[r, niustart + p] <- 1
            }
        }
    }

    r <- r + 1
    col[r] <- k
    Q[r, ONEtoK] <- 1
    ## for (i in 1:k) Q[r, i] <- 1

    for (i in ONEtoN) {
        r <- r + 1
        for (p in ONEtoK) {
            ## d <- st[[p]] # not needed
            ipos <- match(labels[i], ROWNAMES[[p]])
            if (!is.na(ipos)) Q[r, astart[p] + ipos] <- 1
        }
    }

    for (p in 1:(k - 1)) {
        r <- r + 1
        for (i in ONEtoN) {
            ## d <- st[[p]]
            ipos <- match(labels[i], ROWNAMES[[p]])
            if (!is.na(ipos)) Q[r, astart[p] + ipos] <- 1
        }
    }
    a <- solve(Q, col, 1e-19)
    for (i in ONEtoN) {
        for (j in ONEtoN) {
            if (i == j) {
                X[i, j] <- V[i, j] <- 0
                next
            }
            sum <- 0
            sumv <- 0
            for (p in ONEtoK) {
                d <- st[[p]]
                pos <- match(labels[c(i, j)], ROWNAMES[[p]])
                if (all(!is.na(pos))) {
                    ipos <- pos[1L]
                    jpos <- pos[2L]
                    dij <- d[ipos, jpos]
                    sum <- sum + sp[p] * (a[p] * dij + a[astart[p] + ipos] + a[astart[p] + jpos])
                    sumv <- sumv + sp[p] * (a[p] * dij)^2
                }
            }
            X[i, j] <- sum / w[i, j]
            V[i, j] <- sumv / (w[i, j])^2
        }
    }
    list(X, V)
}
