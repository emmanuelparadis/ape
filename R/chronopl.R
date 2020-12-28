## chronopl.R (2012-02-09)

##   Molecular Dating With Penalized Likelihood

## Copyright 2005-2012 Emmanuel Paradis

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

chronopl <-
    function(phy, lambda, age.min = 1, age.max = NULL,
             node = "root", S = 1, tol = 1e-8,
             CV = FALSE, eval.max = 500, iter.max = 500, ...)
{
    n <- length(phy$tip.label)
    ROOT <- n + 1L
    if (identical(node, "root")) node <- ROOT
    if (any(node <= n))
        stop("node numbers should be greater than the number of tips")
    zerobl <- which(phy$edge.length <= 0)
    if (length(zerobl)) {
        if (any(phy$edge[zerobl, 2] <= n))
            stop("at least one terminal branch is of length zero:
  you should remove it to have a meaningful estimation.")
        else {
            warning("at least one internal branch is of length zero:
  it was collapsed and some nodes have been deleted.")
            if (length(node) == 1 && node == ROOT)
                phy <- di2multi(phy)
            else {
                tmp <- FALSE
                if (is.null(phy$node.label)) {
                    tmp <- !tmp
                    phy$node.label <- paste("node", 1:phy$Nnode)
                }
                node.lab <- phy$node.label[node - n]
                phy <- di2multi(phy)
                node <- match(node.lab, phy$node.label) + n
                if (tmp) phy$node.label <- NULL
            }
        }
    }
    m <- phy$Nnode
    el <- phy$edge.length
    e1 <- phy$edge[, 1L]
    e2 <- phy$edge[, 2L]
    N <- length(e1)
    TIPS <- 1:n
    EDGES <- 1:N

    ini.rate <- el
    el <- el/S

    ## `basal' contains the indices of the basal edges
    ## (ie, linked to the root):
    basal <- which(e1 == ROOT)
    Nbasal <- length(basal)

    ## `ind' contains in its 1st column the index of all nonbasal
    ## edges, and in its second column the index of the edges
    ## where these edges come from (ie, this matrix contains pairs
    ## of contiguous edges), eg:

    ##         ___b___    ind:
    ##        |           |   |   |
    ## ___a___|           | b | a |
    ##        |           | c | a |
    ##        |___c___    |   |   |

    ind <- matrix(0L, N - Nbasal, 2)
    ind[, 1] <- EDGES[-basal]
    ind[, 2] <- match(e1[EDGES[-basal]], e2)

    age <- numeric(n + m)

#############################################################################
### This bit sets 'ini.time' and should result in no negative branch lengths

    seq.nod <- .Call("seq_root2tip", phy$edge, n, phy$Nnode, PACKAGE = "ape")

    ini.time <- age
    ini.time[ROOT:(n + m)] <- NA
    ini.time[node] <- if (is.null(age.max)) age.min else (age.min + age.max) / 2

    ## if no age given for the root, find one approximately:
    if (is.na(ini.time[ROOT]))
        ini.time[ROOT] <- if (is.null(age.max)) 3 * max(age.min) else 3 * max(age.max)

    ISnotNA.ALL <- unlist(lapply(seq.nod, function(x) sum(!is.na(ini.time[x]))))
    o <- order(ISnotNA.ALL, decreasing = TRUE)

    for (y in seq.nod[o]) {
        ISNA <- is.na(ini.time[y])
        if (any(ISNA)) {
            i <- 2L # we know the 1st value is not NA, so we start at the 2nd one
            while (i <= length(y)) {
                if (ISNA[i]) { # we stop at the next NA
                    j <- i + 1L
                    while (ISNA[j]) j <- j + 1L # look for the next non-NA
                    nb.val <- j - i
                    by <- (ini.time[y[i - 1L]] - ini.time[y[j]]) / (nb.val + 1)
                    ini.time[y[i:(j - 1L)]] <- ini.time[y[i - 1L]] - by * seq_len(nb.val)
                    i <- j + 1L
                } else i <- i + 1L
            }
        }
    }

    real.edge.length <- ini.time[e1] - ini.time[e2]

    if (any(real.edge.length <= 0))
        stop("some initial branch lengths are zero or negative;
  maybe you need to adjust the given dates -- see '?chronopl' for details")
#############################################################################

    ## because if (!is.null(age.max)), 'node' is modified,
    ## so we copy it in case CV = TRUE:
    node.bak <- node

    ## `unknown.ages' will contain the index of the nodes of unknown age:
    unknown.ages <- n + 1:m

    ## define the bounds for the node ages:
    lower <- rep(tol, length(unknown.ages))
    upper <- rep(1/tol, length(unknown.ages))

    if (!is.null(age.max)) { # are some nodes known within some intervals?
        lower[node - n] <- age.min
        upper[node - n] <- age.max
        ## find nodes known within an interval:
        interv <- which(age.min != age.max)
        ## drop them from the 'node' since they will be estimated:
        node <- node[-interv]
        if (length(node)) age[node] <- age.min[-interv] # update 'age'
    } else age[node] <- age.min

    if (length(node)) {
        unknown.ages <- unknown.ages[n - node] # 'n - node' is simplification for '-(node - n)'
        lower <- lower[n - node]
        upper <- upper[n - node]
    }

    ## `known.ages' contains the index of all nodes (internal and
    ## terminal) of known age:
    known.ages <- c(TIPS, node)

    ## concatenate the bounds for the rates:
    lower <- c(rep(tol, N), lower)
    upper <- c(rep(1 - tol, N), upper)

    minusploglik.gr <- function(rate, node.time) {
        grad <- numeric(N + length(unknown.ages))
        age[unknown.ages] <- node.time
        real.edge.length <- age[e1] - age[e2]
        if (any(real.edge.length < 0)) {
            grad[] <- 0
            return(grad)
        }
        ## gradient for the rates:
        ## the parametric part can be calculated without a loop:
        grad[EDGES] <- real.edge.length - el/rate
        if (Nbasal == 2) { # the simpler formulae if there's a basal dichotomy
            grad[basal[1]] <-
                grad[basal[1]] + lambda*(rate[basal[1]] - rate[basal[2]])
            grad[basal[2]] <-
                grad[basal[2]] + lambda*(rate[basal[2]] - rate[basal[1]])
        } else { # the general case
            for (i in 1:Nbasal)
                grad[basal[i]] <- grad[basal[i]] +
                    lambda*(2*rate[basal[i]]*(1 - 1/Nbasal) -
                            2*sum(rate[basal[-i]])/Nbasal)/(Nbasal - 1)
        }

        for (i in EDGES) {
            ii <- c(which(e2 == e1[i]), which(e1 == e2[i]))
            if (!length(ii)) next
            grad[i] <- grad[i] + lambda*(2*length(ii)*rate[i] - 2*sum(rate[ii]))
        }

        ## gradient for the 'node times'
        for (i in 1:length(unknown.ages)) {
            nd <- unknown.ages[i]
            ii <- which(e1 == nd)
            grad[i + N] <-
                sum(rate[ii] - el[ii]/real.edge.length[ii])#, na.rm = TRUE)
            if (nd != ROOT) {
                ii <- which(e2 == nd)
                grad[i + N] <- grad[i + N] -
                    rate[ii] + el[ii]/real.edge.length[ii]
            }
        }
        grad
    }

    minusploglik <- function(rate, node.time) {
        age[unknown.ages] <- node.time
        real.edge.length <- age[e1] - age[e2]
        if (any(real.edge.length < 0)) return(1e50)
        B <- rate*real.edge.length
        loglik <- sum(-B + el*log(B) - lfactorial(el))
        -(loglik - lambda*(sum((rate[ind[, 1]] - rate[ind[, 2]])^2)
                           + var(rate[basal])))
    }

    out <- nlminb(c(ini.rate, ini.time[unknown.ages]),
                  function(p) minusploglik(p[EDGES], p[-EDGES]),
                  function(p) minusploglik.gr(p[EDGES], p[-EDGES]),
                  control = list(eval.max = eval.max, iter.max = iter.max, ...),
                  lower = lower, upper = upper)

    attr(phy, "ploglik") <- -out$objective
    attr(phy, "rates") <- out$par[EDGES]
    attr(phy, "message") <- out$message
    age[unknown.ages] <- out$par[-EDGES]
    if (CV) ophy <- phy
    phy$edge.length <- age[e1] - age[e2]
    if (CV) attr(phy, "D2") <-
        chronopl.cv(ophy, lambda, age.min, age.max, node.bak,
                    n, S, tol, eval.max, iter.max, ...)
    phy
}

chronopl.cv <- function(ophy, lambda, age.min, age.max, nodes,
                        n, S, tol, eval.max, iter.max, ...)
### ophy: the original phylogeny
### n: number of tips
### Note that we assume here that the order of the nodes
### in node.label are not modified by the drop.tip operation
{
    cat("Doing cross-validation\n")
    BT <- branching.times(ophy)
    D2 <- numeric(n)

    for (i in 1:n) {
        cat("\r  dropping tip ", i, " / ", n, sep = "")
        tr <- drop.tip(ophy, i)
        j <- which(ophy$edge[, 2] == i)
        if (ophy$edge[j, 1] %in% nodes) {
            k <- which(nodes == ophy$edge[j, 1])
            node <- nodes[-k]
            agemin <- age.min[-k]
            agemax <- age.max[-k]
        } else node <- nodes
        if (length(node)) {
            chr <- chronopl(tr, lambda, age.min, age.max, node,
                            S, tol, FALSE, eval.max, iter.max, ...)
            tmp <-
                if (Nnode(chr) == Nnode(ophy)) BT else BT[-(ophy$edge[j, 1] - n)]
            D2[i] <- sum((tmp - branching.times(chr))^2 / tmp)
        } else D2[i] <- 0
    }
    cat("\n")
    D2
}
