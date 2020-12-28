\name{ace}
\alias{ace}
\alias{print.ace}
\alias{logLik.ace}
\alias{deviance.ace}
\alias{AIC.ace}
\alias{anova.ace}
\title{Ancestral Character Estimation}
\description{
  \code{ace} estimates ancestral character states, and the associated
  uncertainty, for continuous and discrete characters. If \code{marginal
    = TRUE}, a marginal estimation procedure is used. With this method,
  the likelihood values at a given node are computed using only the
  information from the tips (and branches) descending from this node.

  The present implementation of marginal reconstruction for discrete
  characters does not calculate the most likely state for each node,
  integrating over all the possible states, over all the other nodes in
  the tree, in proportion to their probability. For more details, see
  the Note below.

  \code{logLik}, \code{deviance}, and \code{AIC} are generic functions
  used to extract the log-likelihood, the deviance, or the Akaike
  information criterion of a fitted object. If no such values are
  available, \code{NULL} is returned.

  \code{anova} is another generic function which is used to compare
  nested models: the significance of the additional parameter(s) is
  tested with likelihood ratio tests. You must ensure that the models
  are effectively nested (if they are not, the results will be
  meaningless). It is better to list the models from the smallest to the
  largest.
}
\usage{
ace(x, phy, type = "continuous", method = if (type == "continuous")
   "REML" else "ML", CI = TRUE,
    model = if (type == "continuous") "BM" else "ER",
    scaled = TRUE, kappa = 1, corStruct = NULL, ip = 0.1,
    use.expm = FALSE, use.eigen = TRUE, marginal = FALSE)
\method{print}{ace}(x, digits = 4, ...)
\method{logLik}{ace}(object, ...)
\method{deviance}{ace}(object, ...)
\method{AIC}{ace}(object, ..., k = 2)
\method{anova}{ace}(object, ...)
}
\arguments{
  \item{x}{a vector or a factor; an object of class \code{"ace"} in the
    case of \code{print}.}
  \item{phy}{an object of class \code{"phylo"}.}
  \item{type}{the variable type; either \code{"continuous"} or
    \code{"discrete"} (or an abbreviation of these).}
  \item{method}{a character specifying the method used for
    estimation. Four choices are possible: \code{"ML"}, \code{"REML"},
    \code{"pic"}, or \code{"GLS"}.}
  \item{CI}{a logical specifying whether to return the 95\% confidence
    intervals of the ancestral state estimates (for continuous
    characters) or the likelihood of the different states (for discrete
    ones).}
  \item{model}{a character specifying the model (ignored if \code{method
      = "GLS"}), or a numeric matrix if \code{type = "discrete"} (see
    details).}
  \item{scaled}{a logical specifying whether to scale the contrast
    estimate (used only if \code{method = "pic"}).}
  \item{kappa}{a positive value giving the exponent transformation of
    the branch lengths (see details).}
  \item{corStruct}{if \code{method = "GLS"}, specifies the correlation
    structure to be used (this also gives the assumed model).}
  \item{ip}{the initial value(s) used for the ML estimation procedure
    when \code{type == "discrete"} (possibly recycled).}
  \item{use.expm}{a logical specifying whether to use the package
    \pkg{expm} to compute the matrix exponential (relevant only if
    \code{type = "d"}). If \code{FALSE}, the function \code{matexpo}
    from \pkg{ape} is used (see details). This option is ignored if
    \code{use.eigen = TRUE} (see next).}
  \item{use.eigen}{a logical (relevant if \code{type = "d"}); if
    \code{TRUE} then the probability matrix is computed with an eigen
    decomposition instead of a matrix exponential (see details).}
  \item{marginal}{a logical (relevant if \code{type = "d"}). By default,
    the joint reconstruction of the ancestral states are done. Set this
    option to \code{TRUE} if you want the marginal reconstruction (see
    details.)}
  \item{digits}{the number of digits to be printed.}
  \item{object}{an object of class \code{"ace"}.}
  \item{k}{a numeric value giving the penalty per estimated parameter;
    the default is \code{k = 2} which is the classical Akaike
    information criterion.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  If \code{type = "continuous"}, the default model is Brownian motion
  where characters evolve randomly following a random walk. This model
  can be fitted by residual maximum likelihood (the default), maximum
  likelihood (Felsenstein 1973, Schluter et al. 1997), least squares
  (\code{method = "pic"}, Felsenstein 1985), or generalized least
  squares (\code{method = "GLS"}, Martins and Hansen 1997, Cunningham et
  al. 1998). In the last case, the specification of \code{phy} and
  \code{model} are actually ignored: it is instead given through a
  correlation structure with the option \code{corStruct}.

  In the setting \code{method = "ML"} and \code{model = "BM"} (this used
  to be the default until \pkg{ape} 3.0-7) the maximum likelihood
  estimation is done simultaneously on the ancestral values and the
  variance of the Brownian motion process; these estimates are then used
  to compute the confidence intervals in the standard way. The REML
  method first estimates the ancestral value at the root (aka, the
  phylogenetic mean), then the variance of the Brownian motion process
  is estimated by optimizing the residual log-likelihood. The ancestral
  values are finally inferred from the likelihood function giving these
  two parameters. If \code{method = "pic"} or \code{"GLS"}, the
  confidence intervals are computed using the expected variances under
  the model, so they depend only on the tree.

  It could be shown that, with a continous character, REML results in
  unbiased estimates of the variance of the Brownian motion process
  while ML gives a downward bias. Therefore the former is recommanded.

  For discrete characters (\code{type = "discrete"}), only maximum
  likelihood estimation is available (Pagel 1994) (see \code{\link{MPR}}
  for an alternative method). The model is specified through a numeric
  matrix with integer values taken as indices of the parameters. The
  numbers of rows and of columns of this matrix must be equal, and are
  taken to give the number of states of the character. For instance,
  \code{matrix(c(0, 1, 1, 0), 2)} will represent a model with two
  character states and equal rates of transition, \code{matrix(c(0, 1,
  2, 0), 2)} a model with unequal rates, \code{matrix(c(0, 1, 1, 1, 0,
  1, 1, 1, 0), 3)} a model with three states and equal rates of
  transition (the diagonal is always ignored). There are short-cuts to
  specify these models: \code{"ER"} is an equal-rates model (e.g., the
  first and third examples above), \code{"ARD"} is an
  all-rates-different model (the second example), and \code{"SYM"} is a
  symmetrical model (e.g., \code{matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0),
  3)}). If a short-cut is used, the number of states is determined from
  the data.

  By default, the likelihood of the different ancestral states of
  discrete characters are computed with a joint estimation procedure
  using a procedure similar to the one described in Pupko et al. (2000).
  If \code{marginal = TRUE}, a marginal estimation procedure is used
  (this was the only choice until ape 3.1-1). With this method, the
  likelihood values at a given node are computed using only the
  information from the tips (and branches) descending from this node.
  With the joint estimation, all information is used for each node. The
  difference between these two methods is further explained in
  Felsenstein (2004, pp. 259-260) and in Yang (2006, pp. 121-126). The
  present implementation of the joint estimation uses a ``two-pass''
  algorithm which is much faster than stochastic mapping while the
  estimates of both methods are very close.

  With discrete characters it is necessary to compute the exponential of
  the rate matrix. The only possibility until \pkg{ape} 3.0-7 was the
  function \code{\link{matexpo}} in \pkg{ape}. If \code{use.expm = TRUE}
  and \code{use.eigen = FALSE}, the function \code{\link[expm]{expm}},
  in the package of the same name, is used. \code{matexpo} is faster but
  quite inaccurate for large and/or asymmetric matrices. In case of
  doubt, use the latter. Since \pkg{ape} 3.0-10, it is possible to use
  an eigen decomposition avoiding the need to compute the matrix
  exponential; see details in Lebl (2013, sect. 3.8.3). This is much
  faster and is now the default.

  Since version 5.2 of \pkg{ape}, \code{ace} can take state uncertainty
  for discrete characters into account: this should be coded with \R's
  \code{\link[base]{NA}} only. More details:

  \url{https://www.mail-archive.com/r-sig-phylo@r-project.org/msg05286.html}
}
\note{
  Liam Revell points out that for discrete characters the ancestral
  likelihood values returned with \code{marginal = FALSE} are actually
  the marginal estimates, while setting \code{marginal = TRUE} returns
  the conditional (scaled) likelihoods of the subtree:

  \url{http://blog.phytools.org/2015/05/about-how-acemarginaltrue-does-not.html}
}
\value{
  an object of class \code{"ace"} with the following elements:

  \item{ace}{if \code{type = "continuous"}, the estimates of the
    ancestral character values.}
  \item{CI95}{if \code{type = "continuous"}, the estimated 95\%
    confidence intervals.}
  \item{sigma2}{if \code{type = "continuous"}, \code{model = "BM"}, and
    \code{method = "ML"}, the maximum likelihood estimate of the
    Brownian parameter.}
  \item{rates}{if \code{type = "discrete"}, the maximum likelihood
    estimates of the transition rates.}
  \item{se}{if \code{type = "discrete"}, the standard-errors of
    estimated rates.}
  \item{index.matrix}{if \code{type = "discrete"}, gives the indices of
    the \code{rates} in the rate matrix.}
  \item{loglik}{if \code{method = "ML"}, the maximum log-likelihood.}
  \item{lik.anc}{if \code{type = "discrete"}, the scaled likelihoods of
    each ancestral state.}
  \item{call}{the function call.}
}
\references{
  Cunningham, C. W., Omland, K. E. and Oakley, T. H. (1998)
  Reconstructing ancestral character states: a critical
  reappraisal. \emph{Trends in Ecology & Evolution}, \bold{13},
  361--366.

  Felsenstein, J. (1973) Maximum likelihood estimation
  of evolutionary trees from continuous characters. \emph{American
  Journal of Human Genetics}, \bold{25}, 471--492.

  Felsenstein, J. (1985) Phylogenies and the comparative
  method. \emph{American Naturalist}, \bold{125}, 1--15.

  Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sunderland:
  Sinauer Associates.

  Lebl, J. (2013) \emph{Notes on Diffy Qs: Differential Equations for
  Engineers}. \url{https://www.jirka.org/diffyqs/}.

  Martins, E. P. and Hansen, T. F. (1997) Phylogenies and the
  comparative method: a general approach to incorporating phylogenetic
  information into the analysis of interspecific data. \emph{American
    Naturalist}, \bold{149}, 646--667.

  Pagel, M. (1994) Detecting correlated evolution on phylogenies: a
  general method for the comparative analysis of discrete
  characters. \emph{Proceedings of the Royal Society of London. Series
    B. Biological Sciences}, \bold{255}, 37--45.

  Pupko, T., Pe'er, I, Shamir, R., and Graur, D. (2000) A fast algorithm
  for joint reconstruction of ancestral amino acid sequences.
  \emph{Molecular Biology and Evolution}, \bold{17}, 890--896.

  Schluter, D., Price, T., Mooers, A. O. and Ludwig, D. (1997)
  Likelihood of ancestor states in adaptive radiation. \emph{Evolution},
  \bold{51}, 1699--1711.

  Yang, Z. (2006) \emph{Computational Molecular Evolution}. Oxford:
  Oxford University Press.
}
\author{Emmanuel Paradis, Ben Bolker}
\seealso{
  \code{\link{MPR}}, \code{\link{corBrownian}}, \code{\link{compar.ou}},
  \code{\link[stats]{anova}}

  Reconstruction of ancestral sequences can be done with the package
  \pkg{phangorn} (see function \code{?ancestral.pml}).
}
\examples{
### Some random data...
data(bird.orders)
x <- rnorm(23)
### Compare the three methods for continuous characters:
ace(x, bird.orders)
ace(x, bird.orders, method = "pic")
ace(x, bird.orders, method = "GLS",
    corStruct = corBrownian(1, bird.orders))
### For discrete characters:
x <- factor(c(rep(0, 5), rep(1, 18)))
ans <- ace(x, bird.orders, type = "d")
#### Showing the likelihoods on each node:
plot(bird.orders, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)
}
\keyword{models}
