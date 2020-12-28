## mcmc.popsize.R (2013-07-19)

##   Run reversible jump MCMC to sample demographic histories

## Copyright 2004-2013 Rainer Opgen-Rhein and Korbinian Strimmer

## Portions of this function are adapted from rjMCMC code by
## Karl W Broman (see http://www.biostat.wisc.edu/~kbroman/)

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

# public function

# run rjMCMC chain
if (getRversion() >= "2.15.1")
    utils::globalVariables(c("loglik", "b.lin", "popsize"))

mcmc.popsize <-
    function(tree, nstep, thinning = 1, burn.in = 0, progress.bar = TRUE,
             method.prior.changepoints = c("hierarchical", "fixed.lambda"),
             max.nodes = 30,
             lambda = 0.5, # "fixed.lambda" method.prior.changepoints
             gamma.shape = 0.5,
             gamma.scale = 2, # gamma distribution from which lambda is drawn (for "hierarchical" method)
             method.prior.heights = c("skyline", "constant", "custom"),
             prior.height.mean, prior.height.var)
{
    method.prior.changepoints <- match.arg(method.prior.changepoints)
    method.prior.heights <- match.arg(method.prior.heights)

    ## Calculate skylineplot, coalescent intervals
    ## and estimated population sizes

    if (inherits(tree, "phylo")) {
        ci <- coalescent.intervals(tree)
        sk1 <- skyline(ci)
    } else if (class(tree) == "coalescentIntervals") {
        ci <- tree
        sk1 <- skyline(ci)
    } else
        stop("tree must be an object of class phylo or coalescentIntervals")

    ## consider possibility of more than one lineage
    ci$lineages <- ci$lineages[sk1$interval.length > 0]
    ci$interval.length <- ci$interval.length[sk1$interval.length > 0]
    data <- sk1$time <- sk1$time[sk1$interval.length > 0]
    sk1$population.size <- sk1$population.size[sk1$interval.length > 0]
    sk1$interval.length <- sk1$interval.length[sk1$interval.length > 0]

    ## constant prior for heights

    if (method.prior.heights == "constant") {
        prior.height.mean <- function(position)
            mean(sk1$population.size)
        prior.height.var <- function(position)
            (mean(sk1$population.size))^2
    }

    ## skyline plot prior for heights

    if (method.prior.heights == "skyline") {
        TIME <- sk1$time
        numb.interv <- 10
        prior.change.times <- abs((0:numb.interv) * max(TIME)/numb.interv)
        prior.height.mean.all <- prior.height.var.all <- vector(length = numb.interv)
        for (p.int in 1:(numb.interv)) {
            left <- p.int
            right <- p.int + 1
            sample.pop <- sk1$population.size[sk1$time >= prior.change.times[left] & sk1$time <= prior.change.times[right]]
            while (length(sample.pop) < 10) {
                if (left > 1) left <- left - 1
                if (right < length(prior.change.times)) right <- right + 1
                sample.pop <- sk1$population.size[sk1$time >= prior.change.times[left] & sk1$time <= prior.change.times[right]]
            }
            prior.height.mean.all[p.int] <- sum(sample.pop)/length(sample.pop)
            prior.height.var.all[p.int] <- sum((sample.pop-prior.height.mean.all[p.int])^2)/(length(sample.pop) - 1)
        }

        prior.height.mean <- function(position) {
            j <- sum(prior.change.times <= position)
            if (j >= length(prior.height.mean.all)) j <- length(prior.height.mean.all)
            prior.mean <- prior.height.mean.all[j]
            prior.mean
        }

        prior.height.var <- function(position) {
            j <- sum(prior.change.times <= position)
            if (j >= length(prior.height.var.all)) j <- length(prior.height.var.all)
            prior.var <- prior.height.var.all[j]
            prior.var
        }
    }

    if (method.prior.heights == "custom") {
        if (missing(prior.height.mean) || missing(prior.height.var))
            stop("custom priors not specified")
    }

    ## set prior
    prior <- vector(length = 4)
    prior[4] <- max.nodes

    ## set initial position of markov chain and likelihood
    pos <- c(0, max(data))
    h <- c(rep(mean(sk1$population.size), 2))

    b.lin <- choose(ci$lineages, 2)
    ## loglik <<- loglik.pop # modified by EP

    ## set lists for data
    count.it <- floor((nstep - burn.in)/thinning)
    save.pos <- save.h <- vector("list", count.it)
    save.loglik <- 1:count.it
    save.steptype <- 1:count.it
    save.accept <- 1:count.it

    ## calculate jump probabilities for given lambda of the prior
    if (method.prior.changepoints == "fixed.lambda") {
        prior[1] <- lambda
        jump.prob <- matrix(ncol = 4, nrow = prior[4] + 1)
        p <- dpois(0:prior[4], prior[1])/ppois(prior[4] + 1, prior[1])
        bk <- c(p[-1]/p[-length(p)], 0)
        bk[bk > 1] <- 1
        dk <- c(0, p[-length(p)]/p[-1])
        dk[dk > 1] <- 1
        mx <- max(bk + dk)
        bk <- bk/mx*0.9
        dk <- dk/mx*0.9
        bk[is.na(bk)] <- 0     # added
        dk[is.na(dk)] <- 0     # added
        jump.prob[, 3] <- bk
        jump.prob[, 4] <- dk
        jump.prob[1, 2] <- 0
        jump.prob[1, 1] <- 1 - bk[1] - dk[1]
        jump.prob[-1, 1] <- jump.prob[-1, 2] <-
            (1 - jump.prob[-1, 3] - jump.prob[-1, 4])/2
    }


    ## calculate starting loglik
    curloglik <- loglik.pop(data, pos, h, b.lin, sk1, ci)

    count.i <- 1

    ## set progress bar
    if (progress.bar == TRUE) {
        dev.new(width = 3, height = 0.7)
        par(mar = c(0.5, 0.5, 2, 0.5))
        plot(x = c(0, 0), y = c(0, 1), type = "l", xlim = c(0, 1), ylim = c(0, 1),
             main = "rjMCMC in progress", ylab = "", xlab = "", xaxs = "i",
             yaxs = "i", xaxt = "n", yaxt = "n")
    }

    ## BEGIN CALCULATION

    for (i in (1:nstep + 1)) {
        if (progress.bar == TRUE) {
            if (i %% 100 == 0) {
                z <- i/nstep
                zt <- (i - 100)/(nstep)
                polygon(c(zt, zt, z, z), c(1, 0, 0, 1), col = "black")
            }
        }

        ## calculate jump probabilities without given lamda
        if (method.prior.changepoints == "hierarchical") {
            prior[1] <- rgamma(1, shape = gamma.shape, scale = gamma.scale)
            jump.prob <- matrix(ncol = 4, nrow = prior[4] + 1)
            p <- dpois(0:prior[4], prior[1]) / ppois(prior[4] + 1, prior[1])
            bk <- c(p[-1]/p[-length(p)], 0)
            bk[bk > 1] <- 1
            dk <- c(0, p[-length(p)]/p[-1])
            dk[dk > 1] <- 1
            mx <- max(bk + dk)
            bk <- bk/mx*0.9
            dk <- dk/mx*0.9
            bk[is.na(bk)] <- 0   # added
            dk[is.na(dk)] <- 0   # added
            jump.prob[, 3] <- bk
            jump.prob[, 4] <- dk
            jump.prob[1, 2] <- 0
            jump.prob[1, 1] <- 1 - bk[1] - dk[1]
            jump.prob[-1, 1] <- jump.prob[-1, 2] <-
                (1 - jump.prob[-1, 3] - jump.prob[-1, 4])/2
        }

        ## determine what type of jump to make
        wh <- sample(1:4, 1, prob = jump.prob[length(h)-1, ])

        if (i %% thinning == 0 & i > burn.in) save.steptype[[count.i]] <- wh

        if (wh == 1) {
            step <- ht.move(data, pos, h, curloglik, prior, b.lin, sk1, ci, prior.height.mean, prior.height.var)
            h <- step[[1]]
            curloglik <- step[[2]]
            if (i %% thinning == 0 & i > burn.in) {
                save.pos[[count.i]] <- pos
                save.h[[count.i]] <- h
                save.loglik[[count.i]] <- step[[2]]
                save.accept[[count.i]] <- step[[3]]
            }
        } else if (wh == 2) {
            step <- pos.move(data, pos, h, curloglik, b.lin, sk1, ci)
            pos <- step[[1]]
            curloglik <- step[[2]]
            if (i %% thinning == 0 & i > burn.in) {
                save.pos[[count.i]] <- pos
                save.h[[count.i]] <- h
                save.loglik[[count.i]] <- step[[2]]
                save.accept[[count.i]] <- step[[3]]
            }
        } else if (wh == 3) {
            step <- birth.step(data, pos, h, curloglik, prior, jump.prob, b.lin, sk1, ci, prior.height.mean, prior.height.var)
            pos <- step[[1]]
            h <- step[[2]]
            curloglik <- step[[3]]
            if (i %% thinning == 0 & i > burn.in) {
                save.pos[[count.i]] <- pos
                save.h[[count.i]] <- h
                save.loglik[[count.i]] <- step[[3]]
                save.accept[[count.i]] <- step[[4]]
            }
        } else {
            step <- death.step(data, pos, h, curloglik, prior, jump.prob, b.lin, sk1, ci, prior.height.mean, prior.height.var)
            pos <- step[[1]]
            h <- step[[2]]
            curloglik <- step[[3]]
            if (i %% thinning == 0 & i > burn.in) {
                save.pos[[count.i]] <- pos
                save.h[[count.i]] <- h
                save.loglik[[count.i]] <- step[[3]]
                save.accept[[count.i]] <- step[[4]]
            }
        }
        if (i %% thinning == 0 & i > burn.in) count.i <- count.i + 1
    }

    if (progress.bar == TRUE) dev.off()

    list(pos = save.pos, h = save.h, loglik = save.loglik,
         steptype = save.steptype, accept = save.accept)
}

## private functions

ht.move <-
    function(data, pos, h, curloglik, prior, b.lin, sk1, ci, prior.height.mean, prior.height.var)
{
    j <- sample(1:length(h), 1)

    prior.mean <- prior.height.mean(pos[j])
    prior.var <- prior.height.var(pos[j])

    prior[3] <- prior.mean/prior.var
    prior[2] <- (prior.mean^2)/prior.var

    newh <- h
    newh[j] <- h[j] * exp(runif(1, -0.5, 0.5))

    newloglik <- loglik.pop(data, pos, newh, b.lin, sk1, ci)
    lr <- newloglik - curloglik

  ratio <- exp(lr + prior[2] * (log(newh[j]) - log(h[j])) - prior[3] * (newh[j] - h[j]))

    if (runif(1, 0, 1) < ratio)
        return(list(newh, newloglik, 1))
    else
        return(list(h, curloglik, 0))
}

pos.move <- function(data, pos, h, curloglik, b.lin, sk1, ci)
{
    j <- if (length(pos) == 3) 2 else sample(2:(length(pos)-1), 1)
    newpos <- pos
    left <- pos[j - 1]
    right <- pos[j + 1]
    newpos[j] <- runif(1, left, right)

    newloglik <- loglik.pop(data, newpos, h, b.lin, sk1, ci)
    lr <- newloglik - curloglik

    ratio <- exp(lr) * (right - newpos[j])*(newpos[j]- left)/
        (right - pos[j])/(pos[j] - left)

    if (runif(1, 0, 1) < ratio)
        return(list(newpos, newloglik, 1))
    else
        return(list(pos, curloglik, 0))
}

birth.step <-
    function(data, pos, h, curloglik, prior, jump.prob, b.lin, sk1, ci, prior.height.mean, prior.height.var)
{
    newpos <- runif(1, 0, pos[length(pos)])
    j <- sum(pos < newpos)

    left <- pos[j]
    right <- pos[j + 1]

    prior.mean <- prior.height.mean(pos[j])
    prior.var <- prior.height.var(pos[j])
    prior[3] <- prior.mean/prior.var
    prior[2] <- (prior.mean^2)/prior.var

    u <- runif(1, -0.5, 0.5)
    oldh <- (((newpos - left)/(right - left))*(h[j + 1] - h[j]) + h[j])
    newheight <- oldh*(1 + u)

    ## ratio
    ## recall that prior = (lambda, alpha, beta, maxk)
    k <- length(pos) - 2
    L <- max(pos)

    prior.logratio <- log(prior[1]) - log(k+1) +  log((2*k + 3)*(2*k + 2)) - 2*log(L) +
        log(newpos - left) + log(right - newpos) - log(right - left) +
            prior[2]*log(prior[3]) - lgamma(prior[2]) +
                (prior[2] - 1) * log(newheight) + prior[3]*(newheight)

    proposal.ratio <- jump.prob[k + 2, 4]*L/jump.prob[k + 1, 3]/(k + 1)
    jacobian <- (((newpos - left)/(right - left))*(h[j + 1] - h[j])) + h[j]

    ## form new parameters
    newpos <- sort(c(pos, newpos))
    newh <- c(h[1:j], newheight, h[(j + 1):length(h)])

    newloglik <- loglik.pop(data, newpos, newh, b.lin, sk1, ci)
    lr <- newloglik - curloglik
    ratio <- exp(lr + prior.logratio) * proposal.ratio * jacobian

    if (runif(1, 0, 1) < ratio)
        return(list(newpos, newh, newloglik, 1))
    else
        return(list(pos, h, curloglik, 0))
}

death.step <-
    function(data, pos, h, curloglik, prior, jump.prob, b.lin, sk1, ci, prior.height.mean, prior.height.var)
{
    ## position to drop
    if (length(pos) == 3) j <- 2
    else j <- sample(2:(length(pos) - 1), 1)

    left <- pos[j - 1]
    right <- pos[j + 1]

    prior.mean <- prior.height.mean(pos[j])
    prior.var <- prior.height.var(pos[j])
    prior[3] <- prior.mean/prior.var
    prior[2] <- (prior.mean^2)/prior.var

    ## get new height
    h.left <- h[j - 1]
    h.right <- h[j + 1]
    newheight <- (((pos[j] - left)/(right - left))*(h.right - h.left) + h.left)

    ## ratio
    ## recall that prior = (lambda, alpha, beta, maxk)
    k <- length(pos) - 3
    L <- max(pos)

    prior.logratio <- log(k+1) - log(prior[1]) -  log(2*(k + 1)*(2*k + 3)) + 2*log(L) -
        log(pos[j] - left) - log(right - pos[j]) + log(right - left) -
            prior[2]*log(prior[3]) + lgamma(prior[2]) -
                (prior[2]-1) * log(newheight) - prior[3]*(newheight)
    proposal.ratio <- (k + 1)*jump.prob[k + 1, 3]/jump.prob[k + 2, 4]/L
    jacobian <- ((pos[j] - left)/(right - left))*(h[j + 1] - h[j - 1]) + h[j - 1]

    ## form new parameters
    newpos <- pos[-j]
    newh <- h[-j]

    newloglik <- loglik.pop(data, newpos, newh, b.lin, sk1, ci)
    lr <- newloglik - curloglik
    ratio <- exp(lr + prior.logratio) * proposal.ratio * (jacobian^(-1))

    if (runif(1, 0, 1) < ratio)
        return(list(newpos, newh, newloglik, 1))
    else
        return(list(pos, h, curloglik, 0))


}

# calculate the log likelihood for a set of data
loglik.pop <-
    function(time = sk1$time, pos = c(0, max(sk1$time)), h = mean(sk1$population.size), b = b.lin, sk1, ci)
{
    data.time <- c(0, time)

    leftside <- 0
    i <- 1
    h1 <- c(h, h[length(h)])
    pos1 <- c(pos, pos[length(pos)])
    while (i < length(time)) {
        left.pos <- sum(data.time[i + 1] >= pos)
        right.pos <- left.pos + 1
        h.mix <- (((data.time[i + 1] - pos[left.pos])/(pos[right.pos] - pos[left.pos]))*(h[right.pos] - h[left.pos])) + h[left.pos]
        leftside <- leftside + log(b[i]/h.mix)
        i <- i + 1
    }

    rightside <- 0
    time1 <- c(0, time)
    time.count <- 1

    ## heigths of jumps
    jumps <- sort(c(time1, pos))
    h.jumps <- jumps
    while (time.count <= length(jumps)) {
        left.pos <- sum(jumps[time.count] >= pos)
        right.pos <- left.pos + 1
        h.jumps[time.count] <- (((jumps[time.count] - pos[left.pos])/(pos[right.pos] - pos[left.pos]))*(h[right.pos] - h[left.pos])) + h[left.pos]
        if (is.na(h.jumps[time.count])) h.jumps[time.count] <- h[left.pos]
        time.count <- time.count + 1
    }

    ## Vector for lineages
    i <- 1
    lineages.jumps <- jumps
    while (i <= length(jumps)) {
        lineages.jumps[i] <- sum(jumps[i] >= time)
        if (lineages.jumps[i] == 0) lineages.jumps[i] <- 1
        i <- i + 1
    }
    lineage <- ci$lineages[lineages.jumps]
    b1 <- choose(lineage, 2)

    ## Integral
    a <- (h.jumps[-1] - h.jumps[-length(h.jumps)])/(jumps[-1] - jumps[-length(jumps)])
    c <- h.jumps[-1] - jumps[-1] * a
    area <- (1/a) * log(a*jumps[-1] + c) - (1/a)*log(a * jumps[-length(jumps)] + c)
    stepfunction <- (jumps[-1] - jumps[-length(jumps)])/h.jumps[-1]
    area[is.na(area)] <- stepfunction[is.na(area)]

    rightside <- sum(area * b1[-1])

    loglik <- leftside - rightside
    loglik
}
