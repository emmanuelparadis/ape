\name{binaryPGLMM}
\alias{binaryPGLMM}
\alias{binaryPGLMM.sim}
\alias{print.binaryPGLMM}
\title{Phylogenetic Generalized Linear Mixed Model for Binary Data}
\description{
  binaryPGLMM performs linear regression for binary phylogenetic data, estimating regression coefficients with approximate standard errors. It simultaneously estimates the strength of phylogenetic signal in the residuals and gives an approximate conditional likelihood ratio test for the hypothesis that there is no signal. Therefore, when applied without predictor (independent) variables, it gives a test for phylogenetic signal for binary data. The method uses a GLMM approach, alternating between penalized quasi-likelihood (PQL) to estimate the "mean components" and restricted maximum likelihood (REML) to estimate the "variance components" of the model.

  binaryPGLMM.sim is a companion function that simulates binary phylogenetic data of the same structure analyzed by binaryPGLMM.
}
\usage{
binaryPGLMM(formula, data = list(), phy, s2.init = 0.1,
            B.init = NULL, tol.pql = 10^-6, maxit.pql = 200,
            maxit.reml = 100)

binaryPGLMM.sim(formula, data = list(), phy, s2 = NULL, B = NULL, nrep = 1)

\method{print}{binaryPGLMM}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{formula}{a two-sided linear formula object describing the
    fixed-effects of the model; for example, Y ~ X.}
  \item{data}{a data frame containing the variables named in formula.}
  \item{phy}{a phylogenetic tree as an object of class "phylo".}
  \item{s2.init}{an initial estimate of s2, the scaling component of the
    variance in the PGLMM. A value of s2 = 0 implies no phylogenetic
    signal. Note that the variance-covariance matrix given by the
    phylogeny phy is scaled to have determinant = 1.}
  \item{B.init}{initial estimates of B, the matrix containing regression
    coefficients in the model. This matrix must have
    dim(B.init)=c(p+1,1), where p is the number of predictor
    (independent) variables; the first element of B corresponds to the
    intercept, and the remaining elements correspond in order to the
    predictor (independent) variables in the model.}
  \item{tol.pql}{a control parameter dictating the tolerance for
    convergence for the PQL optimization.}
  \item{maxit.pql}{a control parameter dictating the maximum number of
    iterations for the PQL optimization.}
  \item{maxit.reml}{a control parameter dictating the maximum number of
    iterations for the REML optimization.}
  \item{x}{an object of class "binaryPGLMM".}
  \item{s2}{in binaryPGLMM.sim, value of s2. See s2.init.}
  \item{B}{in binaryPGLMM.sim, value of B, the matrix containing regression coefficients in the model. See B.init.}
  \item{nrep}{in binaryPGLMM.sim,  number of compete data sets produced.}
  \item{digits}{the number of digits to print.}
  \item{\dots}{further arguments passed to \code{print}.}
}
\details{
     The function estimates parameters for the model

  \deqn{Pr(Y = 1) = q }
  \deqn{q = inverse.logit(b0 + b1 * x1 + b2 * x2 + \dots + \epsilon)}
  \deqn{\epsilon ~ Gaussian(0, s2 * V) }

where \eqn{V} is a variance-covariance matrix derived from a phylogeny (typically under the assumption of Brownian motion evolution). Although mathematically there is no requirement for \eqn{V} to be ultrametric, forcing \eqn{V} into ultrametric form can aide in the interpretation of the model, because in regression for binary dependent variables, only the off-diagonal elements (i.e., covariances) of matrix \eqn{V} are biologically meaningful (see Ives & Garland 2014).

The function converts a phylo tree object into a variance-covariance matrix, and further standardizes this matrix to have determinant = 1. This in effect standardizes the interpretation of the scalar s2. Although mathematically not required, it is a very good idea to standardize the predictor (independent) variables to have mean 0 and variance 1. This will make the function more robust and improve the interpretation of the regression coefficients. For categorical (factor) predictor variables, you will need to construct 0-1 dummy variables, and these should not be standardized (for obvious reasons).

The estimation method alternates between PQL to obtain estimates of the mean components of the model (this is the standard approach to estimating GLMs) and REML to obtain estimates of the variance components. This method gives relatively fast and robust estimation. Nonetheless, the estimates of the coefficients B will generally be upwards bias, as is typical of estimation for binary data. The standard errors of B are computed from the PQL results conditional on the estimate of s2 and therefore should tend to be too small. The function returns an approximate P-value for the hypothesis of no phylogenetic signal in the residuals (i.e., H0:s2 = 0) using an approximate likelihood ratio test based on the conditional REML likelihood (rather than the marginal likelihood). Simulations have shown that these P-values tend to be high (giving type II errors: failing to identify variances that in fact are statistically significantly different from zero).

It is a good idea to confirm statistical inferences using parametric bootstrapping, and the companion function binaryPGLMM.sim gives a simply tool for this. See Examples below.
}
\value{
  An object of class "binaryPGLMM".

  \item{formula}{formula specifying the regression model.}
  \item{B}{estimates of the regression coefficients.}
  \item{B.se}{approximate PQL standard errors of the regression
    coefficients.}
  \item{B.cov}{approximate PQL covariance matrix for the regression
    coefficients.}
  \item{B.zscore}{approximate PQL Z scores for the regression
    coefficients.}
  \item{B.pvalue}{approximate PQL tests for the regression coefficients
    being different from zero.}
  \item{s2}{phylogenetic signal measured as the scalar magnitude of the
    phylogenetic variance-covariance matrix s2 * V.}
  \item{P.H0.s2}{approximate likelihood ratio test of the hypothesis H0
    that s2 = 0. This test is based on the conditional REML (keeping the
    regression coefficients fixed) and is prone to inflated type 1 errors.}
  \item{mu}{for each data point y, the estimate of p that y = 1.}
  \item{b}{for each data point y, the estimate of inverse.logit(p).}
  \item{X}{the predictor (independent) variables returned in matrix form
    (including 1s in the first column).}
  \item{H}{residuals of the form b + (Y - mu)/(mu * (1 - mu)).}
  \item{B.init}{the user-provided initial estimates of B. If B.init is
    not provided, these are estimated using glm() assuming no phylogenetic
    signal. The glm() estimates can generate convergence problems, so
    using small values (e.g., 0.01) is more robust but slower.}
  \item{VCV}{the standardized phylogenetic variance-covariance matrix.}
  \item{V}{estimate of the covariance matrix of H.}
  \item{convergeflag}{flag for cases when convergence failed.}
  \item{iteration}{number of total iterations performed.}
  \item{converge.test.B}{final tolerance for B.}
  \item{converge.test.s2}{final tolerance for s2.}
  \item{rcondflag}{number of times B is reset to 0.01. This is done when
    rcond(V) < 10^(-10), which implies that V cannot be inverted.}
  \item{Y}{in binaryPGLMM.sim, the simulated values of Y.}
}
\author{Anthony R. Ives}
\references{
  Ives, A. R. and Helmus, M. R. (2011) Generalized linear mixed models
  for phylogenetic analyses of community structure. \emph{Ecological
    Monographs}, \bold{81}, 511--525.

  Ives, A. R. and Garland, T., Jr. (2014) Phylogenetic regression for
  binary dependent variables. Pages 231--261 \emph{in} L. Z. Garamszegi,
  editor. \emph{Modern Phylogenetic Comparative Methods and Their
  Application in Evolutionary Biology}. Springer-Verlag, Berlin
  Heidelberg.

}
\seealso{
  package \pkg{pez} and its function \code{communityPGLMM};
  package \pkg{phylolm} and its function \code{phyloglm};
  package \pkg{MCMCglmm}
}
\examples{
## Illustration of binaryPGLMM() with simulated data

# Generate random phylogeny

n <- 100
phy <- compute.brlen(rtree(n=n), method = "Grafen", power = 1)

# Generate random data and standardize to have mean 0 and variance 1
X1 <- rTraitCont(phy, model = "BM", sigma = 1)
X1 <- (X1 - mean(X1))/var(X1)

# Simulate binary Y
sim.dat <- data.frame(Y=array(0, dim=n), X1=X1, row.names=phy$tip.label)
sim.dat$Y <- binaryPGLMM.sim(Y ~ X1, phy=phy, data=sim.dat, s2=.5,
                             B=matrix(c(0,.25),nrow=2,ncol=1), nrep=1)$Y

# Fit model
binaryPGLMM(Y ~ X1, phy=phy, data=sim.dat)

\dontrun{
# Compare with phyloglm()
library(phylolm)
summary(phyloglm(Y ~ X1, phy=phy, data=sim.dat))

# Compare with glm() that does not account for phylogeny
summary(glm(Y ~ X1, data=sim.dat, family="binomial"))

# Compare with logistf() that does not account
# for phylogeny but is less biased than glm()
library(logistf)
logistf(Y ~ X1, data=sim.dat)

# Compare with MCMCglmm
library(MCMCglmm)

V <- vcv(phy)
V <- V/max(V)
detV <- exp(determinant(V)$modulus[1])
V <- V/detV^(1/n)

invV <- Matrix(solve(V),sparse=T)
sim.dat$species <- phy$tip.label
rownames(invV) <- sim.dat$species

nitt <- 43000
thin <- 10
burnin <- 3000

prior <- list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)))
summary(MCMCglmm(Y ~ X1, random=~species, ginvers=list(species=invV),
    data=sim.dat, slice=TRUE, nitt=nitt, thin=thin, burnin=burnin,
    family="categorical", prior=prior, verbose=FALSE))

## Examine bias in estimates of B1 and s2 from binaryPGLMM with
# simulated data. Note that this will take a while.

Reps = 1000

s2 <- 0.4
B1 <- 1

meanEsts <- data.frame(n = Inf, B1 = B1, s2 = s2, Pr.s2 = 1, propconverged = 1)

for (n in c(160, 80, 40, 20)) {

  meanEsts.n <- data.frame(B1 = 0, s2 = 0, Pr.s2 = 0, convergefailure = 0)
    for (rep in 1:Reps) {
      phy <- compute.brlen(rtree(n = n), method = "Grafen", power = 1)
      X <- rTraitCont(phy, model = "BM", sigma = 1)
      X <- (X - mean(X))/var(X)

      sim.dat <- data.frame(Y = array(0, dim = n), X = X, row.names = phy$tip.label)
      sim <- binaryPGLMM.sim(Y ~ 1 + X, phy = phy, data = sim.dat, s2 = s2,
                                       B = matrix(c(0,B1), nrow = 2, ncol = 1), nrep = 1)
      sim.dat$Y <- sim$Y

      z <- binaryPGLMM(Y ~ 1 + X, phy = phy, data = sim.dat)

      meanEsts.n[rep, ] <- c(z$B[2], z$s2, z$P.H0.s2, z$convergeflag == "converged")
  }
converged <- meanEsts.n[,4]
meanEsts <- rbind(meanEsts,
                  c(n, mean(meanEsts.n[converged==1,1]),
                            mean(meanEsts.n[converged==1,2]),
                            mean(meanEsts.n[converged==1, 3] < 0.05),
                            mean(converged)))
}
meanEsts

# Results output for B1 = 0.5, s2 = 0.4; n-Inf gives the values used to
# simulate the data
#    n       B1        s2      Pr.s2 propconverged
# 1 Inf 1.000000 0.4000000 1.00000000         1.000
# 2 160 1.012719 0.4479946 0.36153072         0.993
# 3  80 1.030876 0.5992027 0.24623116         0.995
# 4  40 1.110201 0.7425203 0.13373860         0.987
# 5  20 1.249886 0.8774708 0.05727377         0.873


## Examine type I errors for estimates of B0 and s2 from binaryPGLMM()
# with simulated data. Note that this will take a while.

Reps = 1000

s2 <- 0
B0 <- 0
B1 <- 0

H0.tests <- data.frame(n = Inf, B0 = B0, s2 = s2, Pr.B0 = .05,
                       Pr.s2 = .05, propconverged = 1)
for (n in c(160, 80, 40, 20)) {

  ests.n <- data.frame(B1 = 0, s2 = 0, Pr.B0 = 0, Pr.s2 = 0, convergefailure = 0)
  for (rep in 1:Reps) {
    phy <- compute.brlen(rtree(n = n), method = "Grafen", power = 1)
    X <- rTraitCont(phy, model = "BM", sigma = 1)
    X <- (X - mean(X))/var(X)

    sim.dat <- data.frame(Y = array(0, dim = n), X = X, row.names = phy$tip.label)
    sim <- binaryPGLMM.sim(Y ~ 1, phy = phy, data = sim.dat, s2 = s2,
                           B = matrix(B0, nrow = 1, ncol = 1), nrep = 1)
    sim.dat$Y <- sim$Y

    z <- binaryPGLMM(Y ~ 1, phy = phy, data = sim.dat)

    ests.n[rep, ] <- c(z$B[1], z$s2, z$B.pvalue, z$P.H0.s2, z$convergeflag == "converged")
  }

converged <- ests.n[,5]
H0.tests <- rbind(H0.tests,
                  c(n, mean(ests.n[converged==1,1]),
                    mean(ests.n[converged==1,2]),
                    mean(ests.n[converged==1, 3] < 0.05),
                    mean(ests.n[converged==1, 4] < 0.05),
                    mean(converged)))
}
H0.tests

# Results for type I errors for B0 = 0 and s2 = 0; n-Inf gives the values
# used to simulate the data. These results show that binaryPGLMM() tends to
# have lower-than-nominal p-values; fewer than 0.05 of the simulated
# data sets have H0:B0=0 and H0:s2=0 rejected at the alpha=0.05 level.
#     n            B0         s2      Pr.B0      Pr.s2 propconverged
# 1 Inf  0.0000000000 0.00000000 0.05000000 0.05000000         1.000
# 2 160 -0.0009350357 0.07273163 0.02802803 0.04804805         0.999
# 3  80 -0.0085831477 0.12205876 0.04004004 0.03403403         0.999
# 4  40  0.0019303847 0.25486307 0.02206620 0.03711133         0.997
# 5  20  0.0181394905 0.45949266 0.02811245 0.03313253         0.996
}}
\keyword{regression}
