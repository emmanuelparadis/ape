\name{lmorigin}
\alias{lmorigin}
\alias{print.lmorigin}
\alias{lmorigin.ex1}
\alias{lmorigin.ex2}
\title{ Multiple regression through the origin }
\description{
Function \code{\link{lmorigin}} computes a multiple linear regression and performs tests of significance of the equation parameters (F-test of R-square and t-tests of regression coefficients) using permutations.

The regression line can be forced through the origin. Testing the significance in that case requires a special permutation procedure. This option was developed for the analysis of independent contrasts, which requires regression through the origin. A permutation test, described by Legendre & Desdevises (2009), is needed to analyze contrasts that are not normally distributed.
}
\usage{
lmorigin(formula, data, origin=TRUE, nperm=999, method=NULL, silent=FALSE)
}

\arguments{
  \item{formula }{ A formula specifying the bivariate model, as in
  \code{\link{lm}} and \code{\link{aov}}. }
  \item{data}{ A data frame containing the two variables specified in the formula. }
  \item{origin}{ \code{origin = TRUE} (default) to compute regression through the origin; \code{origin = FALSE} to compute multiple regression with estimation of the intercept. }
  \item{nperm}{ Number of permutations for the tests. If \code{nperm =
   0}, permutation tests will not be computed. The default value is \code{nperm = 999}. For large data files, the permutation test is rather slow since the permutation procedure is not compiled. }
  \item{method}{ \code{method = "raw"} computes t-tests of the regression coefficients by permutation of the raw data. \code{method = "residuals"} computes t-tests of the regression coefficients by permutation of the residuals of the full model. If \code{method = NULL}, permutation of the raw data is used to test the regression coefficients in regression through the origin; permutation of the residuals of the full model is used to test the regression coefficients in ordinary multiple regression. }
  \item{silent}{ Informative messages and the time to compute the tests will not be written to the \R console if silent=TRUE. Useful when the function is called by a numerical simulation function. }
}

\details{
The permutation F-test of R-square is always done by permutation of the raw data. When there is a single explanatory variable, permutation of the raw data is used for the t-test of the single regression coefficient, whatever the method chosen by the user. The rationale is found in Anderson & Legendre (1999).

The \code{print.lmorigin} function prints out the results of the parametric tests (in all cases) and the results of the permutational tests (when nperm > 0).
}

\value{

  \item{reg }{The regression output object produced by function \code{lm}. }
  \item{p.param.t.2tail }{Parametric probabilities for 2-tailed tests of the regression coefficients. }
  \item{p.param.t.1tail }{Parametric probabilities for 1-tailed tests of the regression coefficients. Each test is carried out in the direction of the sign of the coefficient. }
  \item{p.perm.t.2tail }{Permutational probabilities for 2-tailed tests of the regression coefficients. }
  \item{p.perm.t.1tail }{Permutational probabilities for 1-tailed tests of the regression coefficients. Each test is carried out in the direction of the sign of the coefficient. }
  \item{p.perm.F }{Permutational probability for the F-test of R-square. }
  \item{origin }{TRUE is regression through the origin has been computed, FALSE if multiple regression with estimation of the intercept has been used. }
  \item{nperm }{Number of permutations used in the permutation tests. }
  \item{method }{Permutation method for the t-tests of the regression coefficients: \code{method = "raw"} or \code{method = "residuals"}.  }
  \item{var.names }{Vector containing the names of the variables used in the regression. }
  \item{call }{The function call.}
}

\author{ Pierre Legendre, Universite de Montreal }

\references{
Anderson, M. J. and Legendre, P. (1999) An empirical comparison of permutation methods for tests of partial regression coefficients in a linear model. \emph{Journal of Statistical Computation and Simulation}, \bold{62}, 271--303.

Legendre, P. and Desdevises, Y. (2009) Independent contrasts and regression through the origin. \emph{Journal of Theoretical Biology}, \bold{259}, 727--743.

Sokal, R. R. and Rohlf, F. J. (1995) \emph{Biometry - The principles and
  practice of statistics in biological research. Third edition.} New
  York: W. H. Freeman.
}

\examples{
## Example 1 from Sokal & Rohlf (1995) Table 16.1
## SO2 air pollution in 41 cities of the USA
data(lmorigin.ex1)
out <- lmorigin(SO2 ~ ., data=lmorigin.ex1, origin=FALSE, nperm=99)
out

## Example 2: Contrasts computed on the phylogenetic tree of Lamellodiscus
## parasites. Response variable: non-specificity index (NSI); explanatory
## variable: maximum host size. Data from Table 1 of Legendre & Desdevises
## (2009).
data(lmorigin.ex2)
out <- lmorigin(NSI ~ MaxHostSize, data=lmorigin.ex2, origin=TRUE, nperm=99)
out

## Example 3: random numbers
y <- rnorm(50)
X <- as.data.frame(matrix(rnorm(250),50,5))
out <- lmorigin(y ~ ., data=X, origin=FALSE, nperm=99)
out

}

\keyword{ multivariate }
