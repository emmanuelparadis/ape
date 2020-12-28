## compar.lynch.R (2002-08-28)

##   Lynch's Comparative Method

## Copyright 2002 Julien Claude

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

compar.lynch <- function(x, G, eps = 1e-4)
{
    if (is.vector(x) || is.data.frame(x)) x <- as.matrix(x)
    alea <- runif(1, 0, 1)
    z <- as.vector(x)
    uz <- apply(x, 2, mean)
    vcvz <- var(x)
    vz <- diag(vcvz)
    nsp <- nrow(x)
    k <- ncol(x)
    X1 <- matrix(0, k, k)
    diag(X1) <- 1
    I <- matrix(0, nsp, nsp)
    diag(I) <- 1
    vara <- trvare <- matrix(NA, k, k)
    nsp1 <- rep(1, nsp)
    X <- X1 %x% nsp1
    compteur <- 0
    vara <- A0 <- alea * vcvz
    vare <- E0 <- (1 - alea) * vcvz
    newu <- u0 <- uz
    Ginv <- solve(G)
    V0 <- vcvz %x% G
    a0 <- e0 <- matrix(0, nsp, k)
    a1 <- e1 <- matrix(1, nsp, k)
    while (any(abs((rbind(a1, e1) - rbind(a0, e0))) > eps)) {
        a1 <- a0
        e1 <- e0
	compteur <- compteur + 1
        Rinv <- solve(E0 %x% I)
        Dinv <- solve(A0 %x% G)
        info <- solve(Rinv + Dinv)
        newa <- solve(Rinv + Dinv) %*% Rinv %*% (z - X %*% u0)
        newe <- z - X %*% u0 - newa
        e0 <- mnewe <- matrix(newe, nsp, k)
        a0 <- mnewa <- matrix(newa, nsp, k)

        for (i in 1:k) {
            for (j in 1:k) {
                trvare[i, j] <- sum(diag(info[(((i - 1) * nsp) + 1):(i * nsp),
                                              (((j - 1) * nsp) + 1):(j * nsp)]))}
        }
        vare <- ((nsp - 1) * var(mnewe) + trvare) / nsp

        for (i in 1:k) {
            for (j in 1:k) {
                vara[i, j] <- (t(mnewa[, i]) %*% Ginv %*% mnewa[, j] +
                              sum(diag(Ginv %*%
                                       info[(((i - 1) * nsp) + 1):(i * nsp),
                                            (((j - 1) * nsp) + 1):(j * nsp)]))) / nsp
            }
        }

        newu <- apply(x - mnewa, 2, mean)
	V  <-  vara %x% G + vare %x% I

	p <- (2 * pi)^(-nsp) * det(V)^(-0.5) *
          exp(-0.5 * t(z - (X %*% newu)) %*% solve(V) %*% (z - (X %*% newu)))
        E0 <- vare
        A0 <- vara
        u0 <- newu
    }
    dimnames(vare) <- dimnames(vara)
    list(vare = vare, vara = vara, A = mnewa, E = mnewe, u = newu, lik = log(p))
}
