## yule.time.R (2009-02-20)

##    Fits the Time-Dependent Yule Model

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

yule.time <- function(phy, birth, BIRTH = NULL, root.time = 0,
                      opti = "nlm", start = 0.01)
{
    opti <- pmatch(opti, c("nlm", "nlminb", "optim"))
    if (is.na(opti)) stop("ambiguous argument 'opti'")
    LAMBDA <- function() x
    body(LAMBDA) <- body(birth)
    formals(LAMBDA) <- alist(t=)
    BT <- branching.times(phy)
    T <- BT[1]
    x <- BT[1] - BT + root.time
    m <- phy$Nnode

    paranam <- c(names(formals(birth)))
    np <- length(paranam)
    start <- rep(start, length.out = np)

    ## Foo is always vectorized
    if (is.null(BIRTH)) {
        Foo <- function(x) {
            n <- length(x)
            res <- numeric(n)
            for (i in 1:n)
                res[i] <- integrate(LAMBDA, x[i], T)$value
            res
        }
    } else {
        environment(BIRTH) <- environment()
        Foo <- function(x) BIRTH(T) - BIRTH(x)
    }

    half.dev <- function(p) {
        for (i in 1:np)
            assign(paranam[i], p[i], pos = sys.frame(1))
        root.term <-
            if (is.null(BIRTH)) integrate(LAMBDA, x[1], T)$value
            else BIRTH(T) - BIRTH(x[1])
        sum(Foo(x)) + root.term - sum(log(LAMBDA(x[2:m])))
    }

    switch(opti,
        {
            out <- nlm(half.dev, start, hessian = TRUE)
            est <- out$estimate
            se <- sqrt(diag(solve(out$hessian)))
            loglik <- lfactorial(m) - out$minimum
        },{
            out <- nlminb(start, half.dev)
            est <- out$par
            se <- NULL
            loglik <- lfactorial(m) - out$objective
        },{
            out <- optim(start, half.dev, hessian = TRUE,
                         control = list(maxit = 1000), method = "BFGS")
            est <- out$par
            se <- sqrt(diag(solve(out$hessian)))
            loglik <- lfactorial(m) - out$value
        })
    names(est) <- paranam
    if (!is.null(se)) names(se) <- paranam
    structure(list(estimate = est, se = se, loglik = loglik), class = "yule")
}
