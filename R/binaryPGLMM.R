## binaryPGLMM.R (2015-03-04)

##   Phylogenetic Generalized Linear Mixed Model for Binary Data

## Copyright 2015 Anthony R. Ives

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

binaryPGLMM <- function(formula, data = list(), phy, s2.init = 0.1, B.init = NULL, 
	tol.pql = 10^-6, maxit.pql = 200, maxit.reml = 100) {

	# Begin pglmm.reml
	pglmm.reml <- function(par, tinvW, tH, tVphy, tX) {
		n <- dim(tX)[1]
		p <- dim(tX)[2]

		ss2 <- abs(Re(par))
		Cd <- ss2 * tVphy
		V <- tinvW + Cd

		LL <- 10^10
		if (sum(is.infinite(V)) == 0) { # & rcond(V) < 10^10) {
			if (all(eigen(V)$values > 0)) { #if(rcond(V) > 10^-10 & all(eigen(V)$values > 0)) {
				invV <- solve(V)
				logdetV <- determinant(V)$modulus[1]
				if (is.infinite(logdetV)) {
					cholV <- chol(V)
					logdetV <- 2 * sum(log(diag(chol(V))))
				}
				LL <- logdetV + t(tH) %*% invV %*% tH + determinant(t(tX) %*% 
					invV %*% tX)$modulus[1]
			}
		}
		return(LL)
	}
	# End pglmm.reml
	
	if (!inherits(phy, "phylo")) 
		stop("Object \"phy\" is not of class \"phylo\".")
	if (is.null(phy$edge.length)) 
		stop("The tree has no branch lengths.")
	if (is.null(phy$tip.label)) 
		stop("The tree has no tip labels.")
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)

	mf <- model.frame(formula = formula, data = data)
	if (nrow(mf) != length(phy$tip.label)) 
		stop("Number of rows of the design matrix does not match with length of the tree.")
	if (is.null(rownames(mf))) {
		warning("No tip labels, order assumed to be the same as in the tree.\n")
		data.names = phy$tip.label
	} else data.names = rownames(mf)
	order <- match(data.names, phy$tip.label)
	if (sum(is.na(order)) > 0) {
		warning("Data names do not match with the tip labels.\n")
		rownames(mf) <- data.names
	} else {
		tmp <- mf
		rownames(mf) <- phy$tip.label
		mf[order, ] <- tmp[1:nrow(tmp), ]
	}

	X <- model.matrix(attr(mf, "terms"), data = mf)
	y <- model.response(mf)
	if (sum(!(y %in% c(0, 1)))) {
		stop("PGLMM.binary requires a binary response (dependent variable).")
	}
	if (var(y) == 0) {
		stop("The response (dependent variable) is always 0 or always 1.")
	}

	p <- ncol(X)
	Vphy <- vcv(phy)
	Vphy <- Vphy/max(Vphy)
	Vphy/exp(determinant(Vphy)$modulus[1]/n)

	# Compute initial estimates if not provided assuming no phylogeny
	if (!is.null(B.init) & length(B.init) != p) {
		warning("B.init not correct length, so computed B.init using glm()")
	}
	if (is.null(B.init) | (!is.null(B.init) & length(B.init) != p)) {
		B.init <- t(matrix(glm(formula = formula, data = data, family = "binomial")$coefficients, ncol = p))
	}
	B <- B.init
	s2 <- s2.init
	b <- matrix(0, nrow = n)
	beta <- rbind(B, b)
	mu <- exp(X %*% B)/(1 + exp(X %*% B))

	XX <- cbind(X, diag(1, nrow = n, ncol = n))
	C <- s2 * Vphy

	est.s2 <- s2
	est.B <- B
	oldest.s2 <- 10^6
	oldest.B <- matrix(10^6, nrow = length(est.B))

	iteration <- 0
	exitflag <- 0
	rcondflag <- 0
	while (((t(est.s2 - oldest.s2) %*% (est.s2 - oldest.s2) > tol.pql^2) | 
		(t(est.B - oldest.B) %*% (est.B - oldest.B)/length(B) > tol.pql^2)) & 
		(iteration <= maxit.pql)) {

		iteration <- iteration + 1
		oldest.s2 <- est.s2
		oldest.B <- est.B

		est.B.m <- B
		oldest.B.m <- matrix(10^6, nrow = length(est.B))
		iteration.m <- 0

		# mean component
		while ((t(est.B.m - oldest.B.m) %*% (est.B.m - oldest.B.m)/length(B) > 
			tol.pql^2) & (iteration.m <= maxit.pql)) {

			iteration.m <- iteration.m + 1
			oldest.B.m <- est.B.m
			invW <- diag(as.vector((mu * (1 - mu))^-1))
			V <- invW + C

			# This flags cases in which V has a very high condition number, which will cause solve() to fail.
			if (sum(is.infinite(V)) > 0 | rcond(V) < 10^-10) {
				rcondflag <- rcondflag + 1
				B <- 0 * B.init + 0.001
				b <- matrix(0, nrow = n)
				beta <- rbind(B, b)
				mu <- exp(X %*% B)/(1 + exp(X %*% B))
				oldest.B.m <- matrix(10^6, nrow = length(est.B))
				invW <- diag(as.vector((mu * (1 - mu))^-1))
				V <- invW + C
			}

			invV <- solve(V)
			Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
			denom <- t(X) %*% invV %*% X
			num <- t(X) %*% invV %*% Z
			B <- as.matrix(solve(denom, num))

			b <- C %*% invV %*% (Z - X %*% B)
			beta <- rbind(B, b)
			mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))

			est.B.m <- B
		}

		# variance component
		H <- Z - X %*% B
		opt <- optim(fn = pglmm.reml, par = s2, tinvW = invW, tH = H, tVphy = Vphy, 
			tX = X, method = "BFGS", control = list(factr = 1e+12, maxit = maxit.reml))
		s2 <- abs(opt$par)
		C <- s2 * Vphy

		est.s2 <- s2
		est.B <- B
	}
	convergeflag <- "converged"
	if (iteration >= maxit.pql | rcondflag >= 3) {
		convergeflag <- "Did not converge; try increasing maxit.pql or starting with B.init values of .001"
	}

	converge.test.s2 <- (t(est.s2 - oldest.s2) %*% (est.s2 - oldest.s2))^0.5
	converge.test.B <- (t(est.B - oldest.B) %*% (est.B - oldest.B))^0.5/length(est.B)

	# Extract parameters
	invW <- diag(as.vector((mu * (1 - mu))^-1))
	V <- invW + C
	invV <- solve(V)
	Z <- X %*% B + b + (y - mu)/(mu * (1 - mu))
	denom <- t(X) %*% invV %*% X
	num <- t(X) %*% invV %*% Z
	B <- solve(denom, num)
	b <- C %*% invV %*% (Z - X %*% B)
	beta <- rbind(B, b)
	mu <- exp(XX %*% beta)/(1 + exp(XX %*% beta))
	H <- Z - X %*% B

	B.cov <- solve(t(X) %*% invV %*% X)
	B.se <- as.matrix(diag(B.cov))^0.5
	B.zscore <- B/B.se
	B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)

	LL <- opt$value
	lnlike.cond.reml <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% 
		X)$modulus[1] - 0.5 * LL

	LL0 <- pglmm.reml(par = 0, tinvW = invW, tH = H, tVphy = Vphy, tX = X)
	lnlike.cond.reml0 <- -0.5 * (n - p) * log(2 * pi) + 0.5 * determinant(t(X) %*% 
		X)$modulus[1] - 0.5 * LL0

	P.H0.s2 <- pchisq(2 * (lnlike.cond.reml - lnlike.cond.reml0), df = 1, 
		lower.tail = F)/2

	results <- list(formula = formula, B = B, B.se = B.se, B.cov = B.cov, 
		B.zscore = B.zscore, B.pvalue = B.pvalue, s2 = s2, P.H0.s2 = P.H0.s2, 
		mu = mu, b = b, B.init = B.init, X = X, H = H, VCV = Vphy, V = V, 
		convergeflag = convergeflag, iteration = iteration, converge.test.s2 = converge.test.s2, 
		converge.test.B = converge.test.B, rcondflag = rcondflag)
	class(results) <- "binaryPGLMM"
	results
}

######################################################
######################################################
# binaryPGLMM.sim
######################################################
######################################################

binaryPGLMM.sim <- function(formula, data = list(), phy, s2 = NULL, B = NULL, 
	nrep = 1) {

	if (!inherits(phy, "phylo")) 
		stop("Object \"phy\" is not of class \"phylo\".")
	if (is.null(phy$edge.length)) 
		stop("The tree has no branch lengths.")
	if (is.null(phy$tip.label)) 
		stop("The tree has no tip labels.")
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)

	mf <- model.frame(formula = formula, data = data)
	if (nrow(mf) != length(phy$tip.label)) 
		stop("Number of rows of the design matrix does not match with length of the tree.")
	if (is.null(rownames(mf))) {
		warning("No tip labels, order assumed to be the same as in the tree.\n")
		data.names = phy$tip.label
	} else data.names = rownames(mf)
	order <- match(data.names, phy$tip.label)
	if (sum(is.na(order)) > 0) {
		warning("Data names do not match with the tip labels.\n")
		rownames(mf) <- data.names
	} else {
		tmp <- mf
		rownames(mf) <- phy$tip.label
		mf[order, ] <- tmp[1:nrow(tmp), ]
	}

	if (is.null(s2)) 
		stop("You must specify s2")
	if (is.null(B)) 
		stop("You must specify B")
	X <- model.matrix(attr(mf, "terms"), data = mf)

	n <- nrow(X)
	p <- ncol(X)
	V <- vcv(phy)
	V <- V/max(V)
	V <- vcv(phy)
	V <- V/max(V)
	V/exp(determinant(V)$modulus[1]/n)
	V <- s2 * V

	if (s2 > 10^-8) {
		iD <- t(chol(V))
	} else {
		iD <- matrix(0, nrow = n, ncol = n)
	}
	Y <- matrix(0, nrow = n, ncol = nrep)
	y <- matrix(0, nrow = n, ncol = nrep)
	for (i in 1:nrep) {
		y[, i] <- X %*% B + iD %*% rnorm(n = n)
		p <- 1/(1 + exp(-y[, i]))
		Y[, i] <- as.numeric(runif(n = n) < p)
	}
	results <- list(Y = Y, y = y, X = X, s2 = s2, B = B, V = V)
	return(results)
}

######################################################
######################################################
# print.binaryPGLMM
######################################################
######################################################

print.binaryPGLMM <- function(x, digits = max(3, getOption("digits") - 3), 
	...) {

	cat("\n\nCall:")
	print(x$formula)
	cat("\n")

	cat("Random effect (phylogenetic signal s2):\n")
	w <- data.frame(s2 = x$s2, Pr = x$P.H0.s2)
	print(w, digits = digits)

	cat("\nFixed effects:\n")
	coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, 
		Pvalue = x$B.pvalue)
	printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
	cat("\n")
}
