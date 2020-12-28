\name{corphylo}
\alias{corphylo}
\alias{print.corphylo}
\title{Correlations among Multiple Traits with Phylogenetic Signal}
\description{
  This function calculates Pearson correlation coefficients for multiple continuous traits that may have phylogenetic signal, allowing users to specify measurement error as the standard error of trait values at the tips of the phylogenetic tree.  Phylogenetic signal for each trait is estimated from the data assuming that trait evolution is given by a Ornstein-Uhlenbeck process.  Thus, the function allows the estimation of phylogenetic signal in multiple traits while incorporating correlations among traits. It is also possible to include independent variables (covariates) for each trait to remove possible confounding effects.  corphylo() returns the correlation matrix for trait values, estimates of phylogenetic signal for each trait, and regression coefficients for independent variables affecting each trait.
}
\usage{
corphylo(X, U = list(), SeM = NULL, phy = NULL, REML = TRUE,
method = c("Nelder-Mead", "SANN"), constrain.d = FALSE, reltol = 10^-6,
maxit.NM = 1000, maxit.SA = 1000, temp.SA = 1, tmax.SA = 1, verbose = FALSE)

\method{print}{corphylo}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{X}{a n x p matrix with p columns containing the values for the n taxa. Rows of X should have rownames matching the taxon names in phy.}
  \item{U}{a list of p matrices corresponding to the p columns of X, with each matrix containing independent variables for the corresponding column of X. The rownames of each matrix within U must be the same as X, or alternatively, the order of values in rows must match those in X. If U is omitted, only the mean (aka intercept) for each column of X is estimated. If U[[i]] is NULL, only an intercept is estimated for X[, i]. If all values of U[[i]][j] are the same, this variable is automatically dropped from the analysis (i.e., there is no offset in the regression component of the model).}
  \item{SeM}{a n x p matrix with p columns containing standard errors of the trait values in X. The rownames of SeM must be the same as X, or alternatively, the order of values in rows must match those in X. If SeM is omitted, the trait values are assumed to be known without error. If only some traits have mesurement errors, the remaining traits can be given zero-valued standard errors.}
  \item{phy}{a phylo object giving the phylogenetic tree.  The rownames of phy must be the same as X, or alternatively, the order of values in rows must match those in X.}
  \item{REML}{whether REML or ML is used for model fitting.}
  \item{method}{in optim(), either Nelder-Mead simplex minimization or SANN (simulated annealing) minimization is used. If SANN is used, it is followed by Nelder-Mead minimization.}
  \item{constrain.d}{if constrain.d is TRUE, the estimates of d are constrained to be between zero and 1. This can make estimation more stable and can be tried if convergence is problematic. This does not necessarily lead to loss of generality of the results, because before using corphylo, branch lengths of phy can be transformed so that the "starter" tree has strong phylogenetic signal.}
  \item{reltol}{a control parameter dictating the relative tolerance for convergence in the optimization; see optim().}
  \item{maxit.NM}{a control parameter dictating the maximum number of iterations in the optimization with Nelder-Mead minimization; see optim().}
  \item{maxit.SA}{a control parameter dictating the maximum number of iterations in the optimization with SANN minimization; see optim().}
  \item{temp.SA}{a control parameter dictating the starting temperature in the optimization with SANN minimization; see optim().}
  \item{tmax.SA}{a control parameter dictating the number of function evaluations at each temperature in the optimization with SANN minimization; see optim().}
  \item{verbose}{if TRUE, the model logLik and running estimates of the
    correlation coefficients and values of d are printed each iteration
    during optimization.}
  \item{x}{an objects of class corphylo.}
  \item{digits}{the number of digits to be printed.}
  \item{\dots}{arguments passed to and from other methods.}
}
\details{
     For the case of two variables, the function estimates parameters for the model of the form, for example,

  \deqn{X[1] =  B[1,0] + B[1,1] * u[1,1] + \epsilon[1]}
  \deqn{X[2] =  B[2,0] + B[2,1] * u[2,1] + \epsilon[2]}
  \deqn{\epsilon ~ Gaussian(0, V) }

where \eqn{B[1,0]}, \eqn{B[1,1]}, \eqn{B[2,0]}, and \eqn{B[2,1]} are regression coefficients, and \eqn{V} is a variance-covariance matrix containing the correlation coefficient r, parameters of the OU process \eqn{d1} and \eqn{d2}, and diagonal matrices \eqn{M1} and \eqn{M2} of measurement standard errors for \eqn{X[1]} and \eqn{X[2]}. The matrix \eqn{V} is \eqn{2n x 2n}, with \eqn{n x n} blocks given by

\deqn{V[1,1] = C[1,1](d1) + M1}
\deqn{V[1,2] = C[1,2](d1,d2)}
\deqn{V[2,1] = C[2,1](d1,d2)}
\deqn{V[2,2] = C[2,2](d2) + M2}

where \eqn{C[i,j](d1,d2)} are derived from phy under the assumption of joint OU evolutionary processes for each trait (see Zheng et al. 2009).  This formulation extends in the obvious way to more than two traits.

}
\value{
  An object of class "corphylo".

  \item{cor.matrix}{the p x p matrix of correlation coefficients.}
  \item{d}{values of d from the OU process for each trait.}
  \item{B}{estimates of the regression coefficients, including intercepts. Coefficients are named according to the list U. For example, B1.2 is the coefficient corresponding to U[[1]][, 2], and if column 2 in U[[1]] is named "colname2", then the coefficient will be B1.colname2. Intercepts have the form B1.0.}
  \item{B.se}{standard errors of the regression coefficients.}
  \item{B.cov}{covariance matrix for regression coefficients.}
  \item{B.zscore}{Z scores for the regression coefficients.}
  \item{B.pvalue}{tests for the regression coefficients being different from zero.}
  \item{logLik}{he log likelihood for either the restricted likelihood (REML = TRUE) or the overall likelihood (REML = FALSE).}
  \item{AIC}{AIC for either the restricted likelihood (REML = TRUE) or the overall likelihood (REML = FALSE).}
  \item{BIC}{BIC for either the restricted likelihood (REML = TRUE) or the overall likelihood (REML = FALSE).}
  \item{REML}{whether REML is used rather than ML (TRUE or FALSE).}
  \item{constrain.d}{whether or not values of d were constrained to be between 0 and 1 (TRUE or FALSE).}
  \item{XX}{values of X in vectorized form, with each trait X[, i] standardized to have mean zero and standard deviation one.}
  \item{UU}{design matrix with values in UU corresponding to XX; each variable U[[i]][, j] is standardized to have mean zero and standard deviation one.}
  \item{MM}{vector of measurement standard errors corresponding to XX, with the standard errors suitably standardized.}
  \item{Vphy}{the phylogenetic covariance matrix computed from phy and standardized to have determinant equal to one.}
  \item{R}{covariance matrix of trait values relative to the standardized values of XX.}
  \item{V}{overall estimated covariance matrix of residuals for XX including trait correlations, phylogenetic signal, and measurement error variances. This matrix can be used to simulate data for parametric bootstrapping. See examples.}
  \item{C}{matrix V excluding measurement error variances.}
  \item{convcode}{he convergence code provided by optim().}
  \item{niter}{number of iterations performed by optim().}
}
\author{Anthony R. Ives}
\references{
  Zheng, L., A. R. Ives, T. Garland, B. R. Larget, Y. Yu, and K. F. Cao. 2009. New multivariate tests for phylogenetic signal and trait correlations applied to ecophysiological phenotypes of nine \emph{Manglietia} species. \emph{Functional Ecology} \bold{23}:1059--1069.

}
\examples{
## Simple example using data without correlations or phylogenetic
## signal. This illustrates the structure of the input data.

phy <- rcoal(10, tip.label = 1:10)
X <- matrix(rnorm(20), nrow = 10, ncol = 2)
rownames(X) <- phy$tip.label
U <- list(NULL, matrix(rnorm(10, mean = 10, sd = 4), nrow = 10, ncol = 1))
rownames(U[[2]]) <- phy$tip.label
SeM <- matrix(c(0.2, 0.4), nrow = 10, ncol = 2)
rownames(SeM) <- phy$tip.label

corphylo(X = X, SeM = SeM, U = U, phy = phy, method = "Nelder-Mead")

\dontrun{
## Simulation example for the correlation between two variables. The
## example compares the estimates of the correlation coefficients from
## corphylo when measurement error is incorporated into the analyses with
## three other cases: (i) when measurement error is excluded, (ii) when
## phylogenetic signal is ignored (assuming a "star" phylogeny), and (iii)
## neither measurement error nor phylogenetic signal are included.

## In the simulations, variable 2 is associated with a single
## independent variable. This requires setting up a list U that has 2
## elements: element U[[1]] is NULL and element U[[2]] is a n x 1 vector
## containing simulated values of the independent variable.

# Set up parameter values for simulating data
n <- 50
phy <- rcoal(n, tip.label = 1:n)

R <- matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2)
d <- c(0.3, .95)
B2 <- 1

Se <- c(0.2, 1)
SeM <- matrix(Se, nrow = n, ncol = 2, byrow = T)
rownames(SeM) <- phy$tip.label

# Set up needed matrices for the simulations
p <- length(d)

star <- stree(n)
star$edge.length <- array(1, dim = c(n, 1))
star$tip.label <- phy$tip.label

Vphy <- vcv(phy)
Vphy <- Vphy/max(Vphy)
Vphy <- Vphy/exp(determinant(Vphy)$modulus[1]/n)

tau <- matrix(1, nrow = n, ncol = 1) %*% diag(Vphy) - Vphy
C <- matrix(0, nrow = p * n, ncol = p * n)
for (i in 1:p) for (j in 1:p) {
	Cd <- (d[i]^tau * (d[j]^t(tau)) * (1 - (d[i] * d[j])^Vphy))/(1 - d[i] * d[j])
	C[(n * (i - 1) + 1):(i * n), (n * (j - 1) + 1):(j * n)] <- R[i, j] * Cd
}
MM <- matrix(SeM^2, ncol = 1)
V <- C + diag(as.numeric(MM))

## Perform a Cholesky decomposition of Vphy. This is used to generate
## phylogenetic signal: a vector of independent normal random variables,
## when multiplied by the transpose of the Cholesky deposition of Vphy will
## have covariance matrix equal to Vphy.
iD <- t(chol(V))

# Perform Nrep simulations and collect the results
Nrep <- 100
cor.list <- matrix(0, nrow = Nrep, ncol = 1)
cor.noM.list <- matrix(0, nrow = Nrep, ncol = 1)
cor.noP.list <- matrix(0, nrow = Nrep, ncol = 1)
cor.noMP.list <- matrix(0, nrow = Nrep, ncol = 1)
d.list <- matrix(0, nrow = Nrep, ncol = 2)
d.noM.list <- matrix(0, nrow = Nrep, ncol = 2)
B.list <- matrix(0, nrow = Nrep, ncol = 3)
B.noM.list <- matrix(0, nrow = Nrep, ncol = 3)
B.noP.list <- matrix(0, nrow = Nrep, ncol = 3)
for (rep in 1:Nrep) {
	XX <- iD %*% rnorm(2 * n)
	X <- matrix(XX, nrow = n, ncol = 2)
	rownames(X) <- phy$tip.label

	U <- list(NULL, matrix(rnorm(n, mean = 2, sd = 10), nrow = n, ncol = 1))
	rownames(U[[2]]) <- phy$tip.label
	colnames(U[[2]]) <- "V1"
	X[,2] <- X[,2] + B2[1] * U[[2]][,1] - B2[1] * mean(U[[2]][,1])

	z <- corphylo(X = X, SeM = SeM, U = U, phy = phy, method = "Nelder-Mead")
	z.noM <- corphylo(X = X, U = U, phy = phy, method = "Nelder-Mead")
	z.noP <- corphylo(X = X, SeM = SeM, U = U, phy = star, method = "Nelder-Mead")

	cor.list[rep] <- z$cor.matrix[1, 2]
	cor.noM.list[rep] <- z.noM$cor.matrix[1, 2]
	cor.noP.list[rep] <- z.noP$cor.matrix[1, 2]
	cor.noMP.list[rep] <- cor(cbind(lm(X[,1] ~ 1)$residuals, lm(X[,2] ~ U[[2]])$residuals))[1,2]

	d.list[rep, ] <- z$d
	d.noM.list[rep, ] <- z.noM$d

	B.list[rep, ] <- z$B
	B.noM.list[rep, ] <- z.noM$B
	B.noP.list[rep, ] <- z.noP$B

	show(c(rep, z$convcode, z$cor.matrix[1, 2], z$d))
}
correlation <- rbind(R[1, 2], mean(cor.list), mean(cor.noM.list),
                     mean(cor.noP.list), mean(cor.noMP.list))
rownames(correlation) <- c("True", "With SeM and Phy", "Without SeM",
                           "Without Phy", "Without Phy or SeM")
correlation

signal.d <- rbind(d, colMeans(d.list), colMeans(d.noM.list))
rownames(signal.d) <- c("True", "With SeM and Phy", "Without SeM")
signal.d

est.B <- rbind(c(0, 0, B2), colMeans(B.list), colMeans(B.noM.list),
               colMeans(B.noP.list))
rownames(est.B) <- c("True", "With SeM and Phy", "Without SeM", "Without Phy")
colnames(est.B) <- rownames(z$B)
est.B

# Example simulation output
# correlation
                        # [,1]
# True               0.7000000
# With SeM and Phy   0.7055958
# Without SeM        0.3125253
# Without Phy        0.4054043
# Without Phy or SeM 0.3476589

# signal.d
                     # [,1]      [,2]
# True             0.300000 0.9500000
# With SeM and Phy 0.301513 0.9276663
# Without SeM      0.241319 0.4872675

# est.B
                        # B1.0      B2.0     B2.V1
# True              0.00000000 0.0000000 1.0000000
# With SeM and Phy -0.01285834 0.2807215 0.9963163
# Without SeM       0.01406953 0.3059110 0.9977796
# Without Phy       0.02139281 0.3165731 0.9942140
}}
\keyword{regression}
