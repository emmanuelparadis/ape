## CDF.birth.death.R (2019-11-07)

## Functions to Simulate and Fit Time-Dependent Birth-Death Models

## Copyright 2010-2019 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

integrateTrapeze <- function(FUN, from, to, nint = 10)
## compute an integral with a simple trapeze method
## (apparently, Vectorize doesn't give faster calculation)
{
    x <- seq(from = from, to = to, length.out = nint + 1)
    ## reorganized to minimize the calls to FUN:
    out <- FUN(x[1]) + FUN(x[length(x)])
    for (i in 2:nint) out <- out + 2 * FUN(x[i])
    (x[2] - x[1]) * out/2 # (x[2] - x[1]) is the width of the trapezes
}

## case:
## 1: birth and death rates constant
## 2: no primitive available
## 3: primitives are available
## 4: death rate constant, no primitive available
## 5: birth rate constant, no primitive available
## 6: death rate constant, primitive available for birth(t)
## 7: birth rate constant, primitive available for death(t)

.getCase <- function(birth, death, BIRTH = NULL, DEATH = NULL)
{
    if (is.numeric(birth)) {
        if (is.numeric(death)) 1 else {
            if (is.null(DEATH)) 5 else 7
        }
    } else {
        if (is.numeric(death)) {
            if (is.null(BIRTH)) 4 else 6
        } else if (is.null(BIRTH) || is.null(DEATH)) 2 else 3
    }
}

## if (getRversion() >= "2.15.1") -- R 3.2.0 is required for ape
utils::globalVariables("Tmax")

.getRHOetINT <- function(birth, death, BIRTH = NULL, DEATH = NULL, case, fast)
{
    ## build the RHO(), \rho(t), and INT(), I(t), functions
    switch (case,
        { # case 1:
            RHO <- function(t1, t2) (t2 - t1)*(death - birth)
            INT <- function(t) {
                rho <- death - birth
                death*(exp(rho*(Tmax - t)) - 1)/rho
            }
        },{ # case 2:
            if (fast) {
                RHO <- function(t1, t2)
                    integrateTrapeze(function(t) death(t) - birth(t), t1, t2)
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death(u)
                    integrateTrapeze(FOO, t, Tmax)
                }
            } else {
                RHO <- function(t1, t2)
                    integrate(function(t) death(t) - birth(t), t1, t2)$value
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death(u)
                    integrate(Vectorize(FOO), t, Tmax)$value # Vectorize required
                }
            }
        },{ # case 3:
            RHO <- function(t1, t2)
                DEATH(t2) - BIRTH(t2) - DEATH(t1) + BIRTH(t1)
            INT <- function(t) { # vectorized
                FOO <- function(u) exp(RHO(tt, u)) * death(u)
                out <- t
                for (i in 1:length(t)) {
                    tt <- t[i]
                    out[i] <- integrate(FOO, tt, Tmax)$value
                }
                out
            }
        },{ # case 4:
            if (fast) {
                RHO <- function(t1, t2)
                    death * (t2 - t1) - integrateTrapeze(birth, t1, t2)
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death
                    integrateTrapeze(Vectorize(FOO), t, Tmax)
                }
            } else {
                RHO <- function(t1, t2)
                    death * (t2 - t1) - integrate(birth, t1, t2)$value
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death
                    integrate(Vectorize(FOO), t, Tmax)$value
                }
            }
        },{ # case 5:
            RHO <- function(t1, t2)
                integrate(death, t1, t2)$value - birth * (t2 - t1)
            if (fast) {
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death(u)
                    integrateTrapeze(FOO, t, Tmax)
                }
            } else {
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death(u)
                    integrate(Vectorize(FOO), t, Tmax)$value
                }
            }
        },{ # case 6:
            RHO <- function(t1, t2) death * (t2 - t1) - BIRTH(t2) + BIRTH(t1)
            INT <- function(t) { # vectorized
                FOO <- function(u) exp(RHO(tt, u)) * death
                out <- t
                for (i in 1:length(t)) {
                    tt <- t[i]
                    out[i] <- integrate(FOO, tt, Tmax)$value
                }
                out
            }
        },{ # case 7:
            RHO <- function(t1, t2) DEATH(t2) - DEATH(t1) - birth * (t2 - t1)
            if (fast) {
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death(u)
                    integrateTrapeze(FOO, t, Tmax)
                }
            } else {
                INT <- function(t) {
                    FOO <- function(u) exp(RHO(t, u)) * death(u)
                    integrate(Vectorize(FOO), t, Tmax)$value
                }
            }
        })
    list(RHO, INT)
}

CDF.birth.death <-
    function(birth, death, BIRTH = NULL, DEATH = NULL, Tmax, x, case, fast = FALSE)
{
    ff <- .getRHOetINT(birth, death, BIRTH, DEATH, case, fast)
    RHO <- ff[[1]]
    INT <- ff[[2]]
    environment(INT) <- environment() # so that INT() can find Tmax
    .CDF.birth.death2(RHO, INT, birth, death, BIRTH, DEATH,
                      Tmax, x, case, fast)
}

.CDF.birth.death2 <-
    function(RHO, INT, birth, death, BIRTH, DEATH, Tmax, x, case, fast)
{
    Pi <- if (case %in% c(1, 5, 7))
        function(t) (1/(1 + INT(t)))^2 * 2 * exp(-RHO(0, t)) * birth
    else
        function(t) (1/(1 + INT(t)))^2 * 2 * exp(-RHO(0, t)) * birth(t)

    if (!case %in% c(1, 3, 6)) Pi <- Vectorize(Pi)

    denom <-
        if (fast) integrateTrapeze(Pi, 0, Tmax)
        else integrate(Pi, 0, Tmax)$value
    n <- length(x)
    p <- numeric(n)
    if (fast) {
        for (i in 1:n) p[i] <- integrateTrapeze(Pi, 0, x[i])
    } else {
        for (i in 1:n) p[i] <- integrate(Pi, 0, x[i])$value
    }
    p/denom
}

.makePhylo <- function(edge, edge.length, i)
{
    NODES <- edge > 0
    edge[NODES] <- edge[NODES] + i + 1L
    edge[!NODES] <- 1:(i + 1L)

    phy <- list(edge = edge, edge.length = edge.length,
                tip.label = paste("t", 1:(i + 1), sep = ""), Nnode = i)
    class(phy) <- "phylo"
    attr(phy, "order") <- "cladewise"
    phy
}

rlineage <-
    function(birth, death, Tmax = 50, BIRTH = NULL, DEATH = NULL, eps = 1e-6)
{
    case <- .getCase(birth, death, BIRTH, DEATH)

    rTimeToEvent <- function(t)
    {
        ## CDF of the times to event (speciation or extinction):
        switch (case,
            { # case 1:
                Foo <- function(t, x)
                    1 - exp(-(birth + death)*(x - t))
            },{ # case 2:
                Foo <- function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-integrate(function(t) birth(t) + death(t),
                                       t, x)$value)
                }
            },{ # case 3:
                Foo <- function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-(BIRTH(x) - BIRTH(t) + DEATH(x) - DEATH(t)))
                }
            },{ # case 4:
                Foo <- function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-(integrate(function(t) birth(t), t, x)$value
                              + death*(x - t)))
                }

            },{ # case 5:
                Foo <- function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-(birth*(x - t) +
                              integrate(function(t) death(t), t, x)$value))
                }

            },{ # case 6:
                Foo <- function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-(BIRTH(x) - BIRTH(t) + death*(x - t)))
                }

            },{ # case 7:
                Foo <- function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-(birth*(x - t) + DEATH(x) - DEATH(t)))
                }
            })

        ## generate a random time to event by the inverse method:
        P <- runif(1)
        ## in case speciation probability is so low
        ## that time to speciation is infinite:
        if (Foo(t, Tmax) < P) return(Tmax + 1)
        inc <- 10
        x <- t + inc
        while (inc > eps) { # la precision influe sur le temps de calcul
            if (Foo(t, x) > P) {
                x <- x - inc
                inc <- inc/10
            } else x <- x + inc
        }
        x - t
    }

    if (case == 1)
        speORext <- function(t) birth/(birth + death)
    if (case == 2 || case == 3)
        speORext <- function(t) birth(t)/(birth(t) + death(t))
    if (case == 4 || case == 6)
        speORext <- function(t) birth(t)/(birth(t) + death)
    if (case == 5 || case == 7)
        speORext <- function(t) birth/(birth + death(t))

    ## the recursive function implementing algorithm 1
    foo <- function(node) {
        for (k in 0:1) {
            X <- rTimeToEvent(t[node])
            tmp <- t[node] + X
            ## is the event a speciation or an extinction?
            if (tmp >= Tmax) {
                Y <- 0
                tmp <- Tmax
            } else Y <- rbinom(1, size = 1, prob = speORext(tmp))
            j <<- j + 1L
            edge.length[j] <<- tmp - t[node]
            if (Y) {
                i <<- i + 1L
                t[i] <<- tmp
                ## set internal edge:
                edge[j, ] <<- c(node, i)
                foo(i)
            } else
                ## set terminal edge:
                edge[j, ] <<- c(node, 0L)
        }
    }

    edge <- matrix(0L, 1e5, 2)
    edge.length <- numeric(1e5)
    j <- 0L; i <- 1; t <- 0
    foo(1L)
    .makePhylo(edge[1:j, ], edge.length[1:j], i)
}

drop.fossil <- function(phy, tol = 1e-8)
{
    n <- Ntip(phy)
    x <- dist.nodes(phy)[n + 1, ][1:n]
    drop.tip(phy, which(x < max(x) - tol))
}

rbdtree <-
    function(birth, death, Tmax = 50, BIRTH = NULL, DEATH = NULL, eps = 1e-6)
{
    case <- .getCase(birth, death, BIRTH, DEATH)
    ff <- .getRHOetINT(birth, death, BIRTH, DEATH, case, FALSE)
    RHO <- ff[[1]]
    INT <- ff[[2]]
    ## so that RHO() and INT() can find Tmax:
    environment(RHO) <- environment(INT) <- environment()

    rtimetospe <- function(t)
    {
        ## CDF of the times to speciation:
        Foo <- if (case %in% c(1, 5, 7))
            function(t, x) 1 - exp(-birth*(x - t))
        else {
            if (case %in% c(3, 6))
                function(t, x) 1 - exp(-(BIRTH(x) - BIRTH(t)))
            else {
                function(t, x) {
                    if (t == x) return(0)
                    1 - exp(-integrate(birth, t, x)$value)
                }
            }
        }
        ## generate a random time to speciation by the inverse method:
        P <- runif(1)
        ## in case speciation probability is so low
        ## that time to speciation is infinite:
        if (Foo(t, Tmax) < P) return(Tmax + 1)
        inc <- 10
        x <- t + inc
        while (inc > eps) { # la precision influe sur le temps de calcul
            if (Foo(t, x) > P) {
                x <- x - inc
                inc <- inc/10
            } else x <- x + inc
        }
        x - t
    }

    ## the recursive function implementing algorithm 2
    foo <- function(node, start) {
        node <- node # make a local copy
        for (k in 0:1) {
            tau <- start # because tau is changed below
            NoDesc <- TRUE
            X <- rtimetospe(tau)
            while (X < Tmax - tau) {
                tau <- tau + X
                ## does the new lineage survive until Tmax?
                Y <- rbinom(1, size = 1, prob = 1/(1 + INT(tau)))
                if (Y) {
                    i <<- i + 1L
                    t[i] <<- tau
                    ## set internal edge:
                    j <<- j + 1L
                    edge[j, ] <<- c(node, i)
                    edge.length[j] <<- tau - t[node]
                    foo(i, t[i])
                    NoDesc <- FALSE
                    break
                }
                X <- rtimetospe(tau)
            }
            ## set terminal edge:
            if (NoDesc) {
                j <<- j + 1L
                edge[j, 1] <<- node # the 2nd column is already set to 0
                edge.length[j] <<- Tmax - t[node]
            }
        }
    }

    edge <- matrix(0L, 1e5, 2)
    edge.length <- numeric(1e5)
    j <- 0L; i <- 1L; t <- 0
    foo(1L, 0)
    .makePhylo(edge[1:j, ], edge.length[1:j], i)
}

bd.time <- function(phy, birth, death, BIRTH = NULL, DEATH = NULL,
                    ip, lower, upper, fast = FALSE, boot = 0, trace = 0)
{
    guess.bounds <- if (missing(lower)) TRUE else FALSE
    BIG <- 1e10
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

    case <- .getCase(birth, death, BIRTH, DEATH)

    if (is.function(birth)) {
        paranam <- names(formals(birth))
        if (guess.bounds) {
            upper <- rep(BIG, length(paranam))
            lower <- -upper
        }
        formals(birth) <- alist(t=)
        environment(birth) <- environment()
        if (!is.null(BIRTH)) environment(BIRTH) <- environment()
    } else {
        paranam <- "birth"
        if (guess.bounds) {
            upper <- 1
            lower <- 0
        }
    }

    if (is.function(death)) {
        tmp <- names(formals(death))
        np2 <- length(tmp)
        if (guess.bounds) {
            upper <- c(upper, rep(BIG, np2))
            lower <- c(lower, rep(-BIG, np2))
        }
        paranam <- c(paranam, tmp)
        formals(death) <- alist(t=)
        environment(death) <- environment()
        if (!is.null(DEATH)) environment(DEATH) <- environment()
    } else {
        paranam <- c(paranam, "death")
        if (guess.bounds) {
            upper <- c(upper, .1)
            lower <- c(lower, 0)
        }
    }

    np <- length(paranam)

    ff <- .getRHOetINT(birth, death, BIRTH, DEATH, case = case, fast = fast)
    RHO <- ff[[1]]
    INT <- ff[[2]]
    environment(RHO) <- environment(INT) <- environment()

    x <- branching.times(phy)
    n <- length(x)
    Tmax <- x[1]
    x <- Tmax - x # change the time scale so the root is t=0
    x <- sort(x)

    foo <- function(para) {
        for (i in 1:np)
            assign(paranam[i], para[i], pos = sys.frame(1))
        p <- CDF.birth.death(birth, death, BIRTH, DEATH, Tmax = Tmax,
                             x = x, case = case, fast = fast)
        ## w is the probability of the observed tree size (= number of tips)
        w <- PrNt(0, Tmax, Ntip(phy))
        ## p is the expected CDF of branching times
        ## ecdf(x)(x) is the observed CDF
        sum((1:n/n - p)^2)/w # faster than sum((ecdf(x)(x) - p)^2)/w
    }

    if (missing(ip)) ip <- (upper - lower)/2

    out <- nlminb(ip, foo, control = list(trace = trace, eval.max = 500),
                  upper = upper, lower = lower)

    names(out$par) <- paranam
    names(out)[2] <- "SS"

    if (boot) { # nonparametric version
        PAR <- matrix(NA, boot, np)
        i <- 1L
        while (i <= boot) {
            cat("\rDoing bootstrap no.", i, "\n")
            x <- sort(sample(x, replace = TRUE))
            o <- try(nlminb(ip, foo, control = list(trace = 0, eval.max = 500),
                            upper = upper, lower = lower))
            if (class(o) == "list") {
                PAR[i, ] <- o$par
                i <- i + 1L
            }
        }
        out$boot <- PAR
    }
    out
}

LTT <- function(birth = 0.1, death = 0, N = 100, Tmax = 50, PI = 95,
                scaled = TRUE, eps = 0.1, add = FALSE, backward = TRUE,
                ltt.style = list("black", 1, 1),
                pi.style = list("blue", 1, 2), ...)
{
    case <- .getCase(birth, death, NULL, NULL)
    Time <- seq(0, Tmax, eps)
    F <- CDF.birth.death(birth, death, BIRTH = NULL, DEATH = NULL,
                         Tmax = Tmax, x = Time, case = case, fast = TRUE)
    if (PI) {
        i <- (1 - PI/100)/2
        Flow <- qbinom(i, N - 2, F)
        Fup <- qbinom(1 - i, N - 2, F)
        if (scaled) {
            Flow <- Flow/N
            Fup <- Fup/N
        }
    }
    if (!scaled) F <- F * N
    if (backward) Time <- Time - Tmax
    if (add)
        lines(Time, F, "l", col = ltt.style[[1]], lwd = ltt.style[[2]],
              lty = ltt.style[[3]])
    else
        plot(Time, F, "l", col = ltt.style[[1]], lwd = ltt.style[[2]],
             lty = ltt.style[[3]], ylab = "Number of lineages", ...)
    if (PI)
        lines(c(Time, NA, Time), c(Flow, NA, Fup),
              col = pi.style[[1]], lwd = pi.style[[2]], lty = pi.style[[3]])
}

rphylo <-
    function(n, birth, death, BIRTH = NULL, DEATH = NULL,
             T0 = 50, fossils = FALSE, eps = 1e-6)
{
    case <- .getCase(birth, death, BIRTH, DEATH)

    ## Foo(): CDF of the times to event (speciation or extinction)
    ## rTimeToEvent(): generate a random time to event by the inverse method

    switch(case, { # case 1:
        rTimeToEvent <- function(t) t - rexp(1, N * (birth + death)) # much faster than using Foo()
        speORext <- function(t) birth/(birth + death)
        ## Foo <- function(t, x)
        ##     1 - exp(-N*(birth + death)*(x - t))
    },{ # case 2:
        Foo <- function(t, x) {
            if (t == x) return(0)
            1 - exp(-integrate(function(t) birth(t) + death(t),
                               t, x)$value * N)
        }
        speORext <- function(t) birth(t)/(birth(t) + death(t))
    },{ # case 3:
        Foo <- function(t, x) {
            if (t == x) return(0)
            1 - exp(-N*(BIRTH(x) - BIRTH(t) + DEATH(x) - DEATH(t)))
        }
        speORext <- function(t) birth(t)/(birth(t) + death(t))
    },{ # case 4:
        Foo <- function(t, x) {
            if (t == x) return(0)
            1 - exp(-N*(integrate(function(t) birth(t), t, x)$value
                        + death*(x - t)))
        }
        speORext <- function(t) birth(t)/(birth(t) + death)
    },{ # case 5:
        Foo <- function(t, x) {
            if (t == x) return(0)
            1 - exp(-N*(birth*(x - t) +
                            integrate(function(t) death(t), t, x)$value))
        }
        speORext <- function(t) birth/(birth + death(t))
    },{ # case 6:
        Foo <- function(t, x) {
            if (t == x) return(0)
            1 - exp(-N*(BIRTH(x) - BIRTH(t) + death*(x - t)))
        }
        speORext <- function(t) birth(t)/(birth(t) + death)
    },{ # case 7:
        Foo <- function(t, x) {
            if (t == x) return(0)
            1 - exp(-N*(birth*(x - t) + DEATH(x) - DEATH(t)))
        }
        speORext <- function(t) birth/(birth + death(t))
    })

    if (case != 1) {
        rTimeToEvent <- function(t)
        {
            P <- runif(1)
            inc <- 10
            x <- t - inc
            while (inc > eps) {
                if (Foo(x, t) > P) { # fixed by Niko Yasui (2016-07-06) + fixed 2019-11-07
                    x <- x + inc
                    inc <- inc/10
                }
                x <- x - inc
            }
            x
        }
    }

    storage.mode(n) <- "integer"
    N <- n
    t <- T0
    j <- 0L # number of edges already created
    POOL <- seq_len(N) # initial pool (only tips at start)

if (!fossils) {
    Nedge <- 2L * N - 2L
    nextnode <- 2L * N - 1L
    e1 <- integer(Nedge)
    e2 <- integer(Nedge)
    TIME <- numeric(nextnode) # record the times
    TIME[POOL] <- T0 # the times of the n tips are the present time

    while (j < Nedge) {
        X <- rTimeToEvent(t)
        ## is the event a speciation or an extinction?
        Y <- rbinom(1, size = 1, prob = speORext(X))
        if (Y) { # speciation
            i <- sample.int(N, 2)
            fossil <- POOL[i] == 0
            if (any(fossil)) {
                ## we drop the fossil lineage, or the first one if both are fossils
                POOL <- POOL[-i[which(fossil)[1]]]
            } else { # create a node and an edge
                j <- j + 2L
                k <- c(j - 1, j)
                e1[k] <- nextnode
                e2[k] <- POOL[i]
                TIME[nextnode] <- X
                POOL <- c(POOL[-i], nextnode)
                nextnode <- nextnode - 1L
            }
            N <- N - 1L
        } else { # extinction => create a tip, store it in POOL but don't create an edge
            ## fossil lineages are numbered 0 to find them if Y = 1
            N <- N + 1L
            POOL <- c(POOL, 0L)
        }
        t <- X
    }

    Nnode <- n - 1L

} else { # fossils = TRUE
    nextnode <- -1L # nodes are numbered with negatives
    nexttip <- N + 1L # tips are numbered with positives
    e1 <- integer(1e5)
    e2 <- integer(1e5)
    time.tips <- numeric(1e5) # accessed with positive indices
    time.nodes <- numeric(1e5) # accessed with negative indices
    time.tips[POOL] <- T0 # the times of the n living tips are the present time

    while (N > 1) {
        X <- rTimeToEvent(t)
        ## is the event a speciation or an extinction?
        Y <- rbinom(1, size = 1, prob = speORext(X))
        if (Y) { # speciation => create a node
            i <- sample.int(N, 2)
            j <- j + 2L
            k <- c(j - 1, j)
            e1[k] <- nextnode
            e2[k] <- POOL[i]
            time.nodes[-nextnode] <- X
            POOL <- c(POOL[-i], nextnode)
            nextnode <- nextnode - 1L
            N <- N - 1L
        } else { # extinction => create a tip
            N <- N + 1L
            time.tips[nexttip] <- X
            POOL <- c(POOL, nexttip)
            nexttip <- nexttip + 1L
        }
        t <- X
    }

    n <- nexttip - 1L # update n
    Nnode <- n - 1L
    EDGE <- seq_len(j)
    e1 <- e1[EDGE]
    e2 <- e2[EDGE]
    e1 <- e1 + n + Nnode + 1L # e1 has only nodes...
    NODES <- e2 < 0 # ... so this is needed only on e2
    e2[NODES] <- e2[NODES] + n + Nnode + 1L
    ## concatenate the vectors of times after dropping the extra 0's:
    TIME <- c(time.tips[seq_len(n)], rev(time.nodes[seq_len(Nnode)]))
}
    obj <- list(edge = cbind(e1, e2, deparse.level = 0),
                edge.length = TIME[e2] - TIME[e1],
                tip.label = paste0("t", seq_len(n)), Nnode = Nnode)
    class(obj) <- "phylo"
    reorder(obj)
}
