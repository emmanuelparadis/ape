\name{compar.gee}
\alias{compar.gee}
\alias{print.compar.gee}
\alias{drop1.compar.gee}
\alias{predict.compar.gee}
\title{Comparative Analysis with GEEs}
\description{
  \code{compar.gee} performs the comparative analysis using generalized
  estimating equations as described by Paradis and Claude (2002).

  \code{drop1} tests single effects of a fitted model output from
  \code{compar.gee}.

  \code{predict} returns the predicted (fitted) values of the model.
}
\usage{
compar.gee(formula, data = NULL, family = "gaussian", phy, corStruct,
          scale.fix = FALSE, scale.value = 1)
\method{drop1}{compar.gee}(object, scope, quiet = FALSE, ...)
\method{predict}{compar.gee}(object, newdata = NULL, type = c("link", "response"), ...)
}
\arguments{
  \item{formula}{a formula giving the model to be fitted.}
  \item{data}{the name of the data frame where the variables in
    \code{formula} are to be found; by default, the variables are looked
    for in the global environment.}
  \item{family}{a function specifying the distribution assumed for the
    response; by default a Gaussian distribution (with link identity) is
    assumed (see \code{?family} for details on specifying the
    distribution, and on changing the link function).}
  \item{phy}{an object of class \code{"phylo"} (ignored if
    \code{corStruct} is used).}
  \item{corStruct}{a (phylogenetic) correlation structure.}
  \item{scale.fix}{logical, indicates whether the scale parameter should
    be fixed (TRUE) or estimated (FALSE, the default).}
  \item{scale.value}{if \code{scale.fix = TRUE}, gives the value for the
    scale (default: \code{scale.value = 1}).}
  \item{object}{an object of class \code{"compar.gee"} resulting from
    fitting \code{compar.gee}.}
  \item{scope}{<unused>.}
  \item{quiet}{a logical specifying whether to display a warning message
    about eventual ``marginality principle violation''.}
  \item{newdata}{a data frame with column names matching the variables
    in the formula of the fitted object (see
    \code{\link[stats]{predict}} for details).}
  \item{type}{a character string specifying the type of predicted
    values. By default, the linear (link) prediction is returned.}
  \item{\dots}{further arguments to be passed to \code{drop1}.}
}
\details{
  If a data frame is specified for the argument \code{data}, then its
  rownames are matched to the tip labels of \code{phy}. The user must be
  careful here since the function requires that both series of names
  perfectly match, so this operation may fail if there is a typing or
  syntax error. If both series of names do not match, the values in the
  data frame are taken to be in the same order than the tip labels of
  \code{phy}, and a warning message is issued.

  If \code{data = NULL}, then it is assumed that the variables are in
  the same order than the tip labels of \code{phy}.
}
\note{
  The calculation of the phylogenetic degrees of freedom is likely to be
  approximative for non-Brownian correlation structures (this will be
  refined soon).

  The calculation of the quasilikelihood information criterion (QIC)
  needs to be tested.
}
\value{
  \code{compar.gee} returns an object of class \code{"compar.gee"} with
  the following components:
  \item{call}{the function call, including the formula.}
  \item{effect.assign}{a vector of integers assigning the coefficients
    to the effects (used by \code{drop1}).}
  \item{nobs}{the number of observations.}
  \item{QIC}{the quasilikelihood information criterion as defined by Pan
    (2001).}
  \item{coefficients}{the estimated coefficients (or regression parameters).}
  \item{residuals}{the regression residuals.}
  \item{family}{a character string, the distribution assumed for the response.}
  \item{link}{a character string, the link function used for the mean function.}
  \item{scale}{the scale (or dispersion parameter).}
  \item{W}{the variance-covariance matrix of the estimated coefficients.}
  \item{dfP}{the phylogenetic degrees of freedom (see Paradis and Claude
    for details on this).}

  \code{drop1} returns an object of class \code{"\link[stats]{anova}"}.

  \code{predict} returns a vector or a data frame if \code{newdata} is used.
}
\references{
  Pan, W. (2001) Akaike's information criterion in generalized
  estimating equations. \emph{Biometrics}, \bold{57}, 120--125.

  Paradis, E. and Claude J. (2002) Analysis of comparative data using
  generalized estimating equations. \emph{Journal of theoretical
    Biology}, \bold{218}, 175--185.
}

\author{Emmanuel Paradis}

\seealso{
  \code{\link{read.tree}}, \code{\link{pic}},
  \code{\link{compar.lynch}}, \code{\link[stats]{drop1}}
}
\examples{
### The example in Phylip 3.5c (originally from Lynch 1991)
### (the same analysis than in help(pic)...)
tr <- "((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);"
tree.primates <- read.tree(text = tr)
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
### Both regressions... the results are quite close to those obtained
### with pic().
compar.gee(X ~ Y, phy = tree.primates)
compar.gee(Y ~ X, phy = tree.primates)
### Now do the GEE regressions through the origin: the results are quite
### different!
compar.gee(X ~ Y - 1, phy = tree.primates)
compar.gee(Y ~ X - 1, phy = tree.primates)
}
\keyword{regression}
