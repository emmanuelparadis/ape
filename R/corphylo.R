## corphylo.R (2015-05-01)

##   Ancestral Character Estimation

## Copyright 2015 Anthony R. Ives

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

corphylo <- function(X, U = list(), SeM = NULL, phy = NULL, REML = TRUE, method = c("Nelder-Mead", "SANN"), 
	constrain.d = FALSE, reltol = 10^-6, maxit.NM = 1000, maxit.SA = 1000, temp.SA = 1, tmax.SA = 1, verbose = FALSE) {

	# Begin corphylo.LL
	corphylo.LL <- function(par, XX, UU, MM, tau, Vphy, REML, constrain.d, verbose) {

		n <- nrow(X)
		p <- ncol(X)

		L.elements <- par[1:(p + p * (p - 1)/2)]
		L <- matrix(0, nrow = p, ncol = p)
		L[lower.tri(L, diag = T)] <- L.elements
		R <- t(L) %*% L

		if (constrain.d == TRUE) {
			logit.d <- par[(p + p * (p - 1)/2 + 1):length(par)]
			if (max(abs(logit.d)) > 10) 
				return(10^10)
			d <- 1/(1 + exp(-logit.d))
		} else {
			d <- par[(p + p * (p - 1)/2 + 1):length(par)]
			if (max(d) > 10) 
				return(10^10)
		}

		# OU transform
		C <- matrix(0, nrow = p * n, ncol = p * n)
		for (i in 1:p) for (j in 1:p) {
			Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
			C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
		}

		V <- C + diag(as.numeric(MM))
		if (is.nan(rcond(V)) || rcond(V) < 10^-10) 
			return(10^10)
		iV <- solve(V)
		denom <- t(UU) %*% iV %*% UU
		if (is.nan(rcond(denom)) || rcond(denom) < 10^-10) 
			return(10^10)
		num <- t(UU) %*% iV %*% XX
		B <- solve(denom, num)
		B <- as.matrix(B)
		H <- XX - UU %*% B

		logdetV <- -determinant(iV)$modulus[1]
		if (is.infinite(logdetV)) 
			return(10^10)

		if (REML == TRUE) {
			# REML likelihood function		
			LL <- 0.5 * (logdetV + determinant(t(UU) %*% iV %*% UU)$modulus[1] + t(H) %*% iV %*% H)
		} else {
			# ML likelihood function
			LL <- 0.5 * (logdetV + t(H) %*% iV %*% H)
		}

		if (verbose == T) 
			show(c(as.numeric(LL), par))
		return(as.numeric(LL))
	}
	# End corphylo.LL
	
	# Main program
	if (!inherits(phy, "phylo")) 
		stop("Object \"phy\" is not of class \"phylo\".")
	if (is.null(phy$edge.length)) 
		stop("The tree has no branch lengths.")
	if (is.null(phy$tip.label)) 
		stop("The tree has no tip labels.")
	phy <- reorder(phy, "postorder")
	n <- length(phy$tip.label)

	# Input X
	if (dim(X)[1] != n) 
		stop("Number of rows of the data matrix does not match the length of the tree.")
	if (is.null(rownames(X))) {
		warning("No tip labels on X; order assumed to be the same as in the tree.\n")
		data.names = phy$tip.label
	} else data.names = rownames(X)
	order <- match(data.names, phy$tip.label)
	if (sum(is.na(order)) > 0) {
		warning("Data names do not match with the tip labels.\n")
		rownames(X) <- data.names
	} else {
		temp <- X
		rownames(X) <- phy$tip.label
		X[order, ] <- temp[1:nrow(temp), ]
	}
	p <- dim(X)[2]

	# Input SeM
	if (!is.null(SeM)) {
		if (dim(SeM)[1] != n) 
			stop("Number of rows of the SeM matrix does not match the length of the tree.")
		if (is.null(rownames(SeM))) {
			warning("No tip labels on SeM; order assumed to be the same as in the tree.\n")
			data.names = phy$tip.label
		} else data.names = rownames(SeM)
		order <- match(data.names, phy$tip.label)
		if (sum(is.na(order)) > 0) {
			warning("SeM names do not match with the tip labels.\n")
			rownames(SeM) <- data.names
		} else {
			temp <- SeM
			rownames(SeM) <- phy$tip.label
			SeM[order, ] <- temp[1:nrow(temp), ]
		}
	} else {
		SeM <- matrix(0, nrow = n, ncol = p)
	}

	# Input U
	if (length(U) > 0) {
		if (length(U) != p) 
			stop("Number of elements of list U does not match the number of columns in X.")

		for (i in 1:p) {
			if (!is.null(U[[i]])){
				if (dim(U[[i]])[1] != n) 
					stop("Number of rows of an element of U does not match the tree.")
				if (is.null(rownames(U[[i]]))) {
					warning("No tip labels on U; order assumed to be the same as in the tree.\n")
					data.names = phy$tip.label
				} else data.names = rownames(U[[i]])
				order <- match(data.names, phy$tip.label)
				if (sum(is.na(order)) > 0) {
					warning("U names do not match with the tip labels.\n")
					rownames(U[[i]]) <- data.names
				} else {
					temp <- U[[i]]
					rownames(U[[i]]) <- phy$tip.label
					U[[i]][order, ] <- temp[1:nrow(temp), ]
				}
			} else {
				U[[i]] <- matrix(0, nrow=n, ncol=1)
				rownames(U[[i]]) <- phy$tip.label
			}
		}
	}

	# Standardize all variables
	Xs <- X
	for (i in 1:p) Xs[, i] <- (X[, i] - mean(X[, i]))/sd(X[, i])

	if (!is.null(SeM)) {
		SeMs <- SeM
		for (i in 1:p) SeMs[, i] <- SeM[, i]/sd(X[, i])
	}

	if (length(U) > 0) {
		Us <- U
		for (i in 1:p) for (j in 1:ncol(U[[i]])) {
			if (sd(U[[i]][, j]) > 0) {
				Us[[i]][, j] <- (U[[i]][, j] - mean(U[[i]][, j]))/sd(U[[i]][, j])
			} else {
				Us[[i]][, j] <- U[[i]][, j] - mean(U[[i]][, j])
			}
		}
	}

	# Set up matrices
	Vphy <- vcv(phy)
	Vphy <- Vphy/max(Vphy)
	Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)

	XX <- matrix(as.matrix(Xs), ncol = 1)
	MM <- matrix(as.matrix(SeMs^2), ncol = 1)

	UU <- kronecker(diag(p), matrix(1, nrow = n, ncol = 1))
	if (length(U) > 0) {
		zeros <- 0 * (1:p)
		for (i in 1:p) {
			dd <- zeros
			dd[i] <- 1
			u <- kronecker(dd, as.matrix(Us[[i]]))
			for (j in 1:dim(u)[2]) if (sd(u[, j]) > 0) 
				UU <- cbind(UU, u[, j])
		}
	}

	# Compute initial estimates assuming no phylogeny if not provided
	if (length(U) > 0) {
		eps <- matrix(nrow = n, ncol = p)
		for (i in 1:p) {
			if (ncol(U[[i]]) > 0) {
				u <- as.matrix(Us[[i]])
				z <- lm(Xs[, i] ~ u)
				eps[, i] <- resid(z)
			} else {
				eps[, i] <- Xs[, i] - mean(Xs[, i])
			}
		}
		L <- t(chol(cov(eps)))
	} else {
		L <- t(chol(cov(Xs)))
	}
	L.elements <- L[lower.tri(L, diag = T)]
	par <- c(L.elements, array(0.5, dim = c(1, p)))

	tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy

	if (method == "Nelder-Mead") 
		opt <- optim(fn = corphylo.LL, par = par, XX = XX, UU = UU, MM = MM, tau = tau, Vphy = Vphy, REML = REML, verbose = verbose, constrain.d = constrain.d, method = "Nelder-Mead", control = list(maxit = maxit.NM, reltol = reltol))

	if (method == "SANN") {
		opt <- optim(fn = corphylo.LL, par = par, XX = XX, UU = UU, MM = MM, tau = tau, Vphy = Vphy, REML = REML, 
			verbose = verbose, constrain.d = constrain.d, method = "SANN", control = list(maxit = maxit.SA, 
				temp = temp.SA, tmax = tmax.SA, reltol = reltol))
		par <- opt$par
		opt <- optim(fn = corphylo.LL, par = par, XX = XX, UU = UU, MM = MM, tau = tau, Vphy = Vphy, REML = REML, 
			verbose = verbose, constrain.d = constrain.d, method = "Nelder-Mead", control = list(maxit = maxit.NM, 
				reltol = reltol))
	}

	# Extract parameters
	par <- Re(opt$par)
	LL <- opt$value

	L.elements <- par[1:(p + p * (p - 1)/2)]
	L <- matrix(0, nrow = p, ncol = p)
	L[lower.tri(L, diag = T)] <- L.elements
	R <- t(L) %*% L
	Rd <- diag(diag(R)^-0.5)
	cor.matrix <- Rd %*% R %*% Rd

	if (constrain.d == TRUE) {
		logit.d <- par[(p + p * (p - 1)/2 + 1):length(par)]
		d <- 1/(1 + exp(-logit.d))
	} else {
		d <- par[(p + p * (p - 1)/2 + 1):length(par)]
	}

	# OU transform
	C <- matrix(0, nrow = p * n, ncol = p * n)
	for (i in 1:p) for (j in 1:p) {
		Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
		C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
	}

	V <- C + diag(MM)
	iV <- solve(V)
	denom <- t(UU) %*% iV %*% UU
	num <- t(UU) %*% iV %*% XX
	B <- solve(denom, num)
	B <- as.matrix(B)
	B.cov <- solve(t(UU) %*% iV %*% UU)
	H <- XX - UU %*% B

	# Back-transform B
	counter <- 0
	sd.list <- matrix(0, nrow = dim(UU)[2], ncol = 1)
	for (i in 1:p) {
		counter <- counter + 1
		B[counter] <- B[counter] + mean(X[, i])
		sd.list[counter] <- sd(X[, i])
		if (length(U) > 0) {
			for (j in 1:ncol(U[[i]])) {
				if (sd(U[[i]][, j]) > 0) {
					counter <- counter + 1
					B[counter] <- B[counter] * sd(X[, i])/sd(U[[i]][, j])
					sd.list[counter] <- sd(X[, i])/sd(U[[i]][, j])
				}
			}
		}
	}
	B.cov <- diag(as.numeric(sd.list)) %*% B.cov %*% diag(as.numeric(sd.list))
	B.se <- as.matrix(diag(B.cov))^0.5
	B.zscore <- B/B.se
	B.pvalue <- 2 * pnorm(abs(B/B.se), lower.tail = FALSE)

	# RowNames for B
	if (length(U) > 0) {
		B.rownames <- NULL
		for (i in 1:p) {
			B.rownames <- c(B.rownames, paste("B", i, ".0", sep = ""))
			if (ncol(U[[i]]) > 0) 
				for (j in 1:ncol(U[[i]])) if (sd(U[[i]][, j]) > 0) {
					if (is.null(colnames(U[[i]])[j])) 
						B.rownames <- c(B.rownames, paste("B", i, ".", j, sep = ""))
					if (!is.null(colnames(U[[i]])[j])) 
						B.rownames <- c(B.rownames, paste("B", i, ".", colnames(U[[i]])[j], sep = ""))
				}
		}
	} else {
		B.rownames <- NULL
		for (i in 1:p) {
			B.rownames <- c(B.rownames, paste("B", i, ".0", sep = ""))
		}
	}
	rownames(B) <- B.rownames
	rownames(B.cov) <- B.rownames
	colnames(B.cov) <- B.rownames
	rownames(B.se) <- B.rownames
	rownames(B.zscore) <- B.rownames
	rownames(B.pvalue) <- B.rownames

	if (REML == TRUE) {
		logLik <- -0.5 * ((n * p) - ncol(UU)) * log(2 * pi) + 0.5 * determinant(t(XX) %*% XX)$modulus[1] - LL
	} else {
		logLik <- -0.5 * (n * p) * log(2 * pi) - LL
	}
	k <- length(par) + ncol(UU)
	AIC <- -2 * logLik + 2 * k
	BIC <- -2 * logLik + k * (log(n) - log(pi))

	results <- list(cor.matrix = cor.matrix, d = d, B = B, B.se = B.se, B.cov = B.cov, B.zscore = B.zscore, 
		B.pvalue = B.pvalue, logLik = logLik, AIC = AIC, BIC = BIC, REML = REML, constrain.d = constrain.d, 
		XX = XX, UU = UU, MM = MM, Vphy = Vphy, R = R, V = V, C = C, convcode = opt$convergence, niter = opt$counts)
	class(results) <- "corphylo"
	return(results)
}

# Printing corphylo objects
print.corphylo <- function(x, digits = max(3, getOption("digits") - 3), ...) {

	cat("Call to corphylo\n\n")

	logLik = x$logLik
	AIC = x$AIC
	BIC = x$BIC

	names(logLik) = "logLik"
	names(AIC) = "AIC"
	names(BIC) = "BIC"
	print(c(logLik, AIC, BIC), digits = digits)

	cat("\ncorrelation matrix:\n")
	rownames(x$cor.matrix) <- 1:dim(x$cor.matrix)[1]
	colnames(x$cor.matrix) <- 1:dim(x$cor.matrix)[1]
	print(x$cor.matrix, digits = digits)

	cat("\nfrom OU process:\n")
	d <- data.frame(d = x$d)
	print(d, digits = digits)
	if (x$constrain.d == TRUE) 
		cat("\nvalues of d constrained to be in [0, 1]\n")

	cat("\ncoefficients:\n")
	coef <- data.frame(Value = x$B, Std.Error = x$B.se, Zscore = x$B.zscore, Pvalue = x$B.pvalue)
	rownames(coef) <- rownames(x$B)
	printCoefmat(coef, P.values = TRUE, has.Pvalue = TRUE)
	cat("\n")

	if (x$convcode != 0) 
		cat("\nWarning: convergence in optim() not reached\n")
}
