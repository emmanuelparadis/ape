\name{yule.cov}
\alias{yule.cov}
\title{Fits the Yule Model With Covariates}
\usage{
yule.cov(phy, formula, data = NULL)
}
\arguments{
  \item{phy}{an object of class \code{"phylo"}.}
  \item{formula}{a formula specifying the model to be fitted.}
  \item{data}{the name of the data frame where the variables in
    \code{formula} are to be found; by default, the variables are looked
    for in the global environment.}
}
\description{
  This function fits by maximum likelihood the Yule model with
  covariates, that is a birth-only model where speciation rate is
  determined by a generalized linear model.
}
\details{
  The model fitted is a generalization of the Yule model where the
  speciation rate is determined by:

  \deqn{\ln\frac{\lambda_i}{1 - \lambda_i} = \beta_1 x_{i1} + \beta_2 x_{i2}
    + \dots + \alpha }{ln(li / (1 - li)) = b1 xi1 + b2 xi2 + ... a}

  where \eqn{\lambda_i}{li} is the speciation rate for species i,
  \eqn{x_{i1}, x_{i2}, \dots}{xi1, xi2, ...} are species-specific
  variables, and \eqn{\beta_1, \beta_2, \dots, \alpha}{b1, b2, ..., a}
  are parameters to be estimated. The term on the left-hand side above
  is a logit function often used in generalized linear models for
  binomial data (see \code{\link[stats]{family}}). The above model can
  be written in matrix form:

  \deqn{\mathrm{logit} \lambda_i = x_i' \beta}{logit li = xi' b}

  The standard-errors of the parameters are computed with the second
  derivatives of the log-likelihood function. (See References for other
  details on the estimation procedure.)

  The function needs three things:

\itemize{
  \item a phylogenetic tree which may contain multichotomies;

  \item a formula which specifies the predictors of the model described
  above: this is given as a standard \R formula and has no response (no
  left-hand side term), for instance: \code{~ x + y}, it can include
  interactions (\code{~ x + a * b}) (see \code{\link[stats]{formula}}
  for details);

  \item the predictors specified in the formula must be accessible to
  the function (either in the global space, or though the \code{data}
  option); they can be numeric vectors or factors. The length and the
  order of these data are important: the number of values (length) must
  be equal to the number of tips of the tree + the number of nodes. The
  order is the following: first the values for the tips in the same
  order than for the labels, then the values for the nodes sequentially
  from the root to the most terminal nodes (i.e., in the order given by
  \code{phy$edge}).
}

The user must obtain the values for the nodes separately.

Note that the method in its present implementation assumes that the
change in a species trait is more or less continuous between two nodes
or between a node and a tip. Thus reconstructing the ancestral values
with a Brownian motion model may be consistent with the present
method. This can be done with the function \code{\link{ace}}.
}
\value{
  A NULL value is returned, the results are simply printed. The output
  includes the deviance of the null (intercept-only) model and a
  likelihood-ratio test of the fitted model against the null model.
  Note that the deviance of the null model is different from the one
  returned by \code{\link{yule}} because of the different parametrizations.
}
\references{
  Paradis, E. (2005) Statistical analysis of diversification with
  species traits. \emph{Evolution}, \bold{59}, 1--12.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{branching.times}}, \code{\link{diversi.gof}},
  \code{\link{diversi.time}}, \code{\link{ltt.plot}},
  \code{\link{birthdeath}}, \code{\link{bd.ext}}, \code{\link{yule}}
}
\examples{
### a simple example with some random data
data(bird.orders)
x <- rnorm(45) # the tree has 23 tips and 22 nodes
### the standard-error for x should be as large as
### the estimated parameter
yule.cov(bird.orders, ~ x)
### another example with a tree that has a multichotomy
data(bird.families)
y <- rnorm(272) # 137 tips + 135 nodes
yule.cov(bird.families, ~ y)
}
\keyword{models}
