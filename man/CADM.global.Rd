\name{CADM.global}
\alias{CADM}
\alias{CADM.global}
\alias{CADM.post}
\title{ Congruence among distance matrices }
\description{
Function \code{\link{CADM.global}} compute and test the coefficient of concordance among several distance matrices through a permutation test.

Function \code{\link{CADM.post}} carries out a posteriori permutation tests of the contributions of individual distance matrices to the overall concordance of the group.

Use in phylogenetic analysis: to identify congruence among distance matrices (D) representing different genes or different types of data. Congruent D matrices correspond to data tables that can be used together in a combined phylogenetic or other type of multivariate analysis.
}
\usage{
CADM.global(Dmat, nmat, n, nperm=99, make.sym=TRUE, weights=NULL,
            silent=FALSE)
CADM.post  (Dmat, nmat, n, nperm=99, make.sym=TRUE, weights=NULL,
             mult="holm", mantel=FALSE, silent=FALSE)
}

\arguments{
  \item{Dmat}{ A text file listing the distance matrices one after the other, with or without blank lines in-between. Each matrix is in the form of a square distance matrix with 0's on the diagonal. }
  \item{nmat}{ Number of distance matrices in file Dmat. }
  \item{n}{ Number of objects in each distance matrix. All matrices must have the same number of objects. }
  \item{nperm}{ Number of permutations for the tests of significance. }
  \item{make.sym}{ TRUE: turn asymmetric matrices into symmetric matrices by averaging the two triangular portions. FALSE: analyse asymmetric matrices as they are. }
  \item{weights}{ A vector of positive weights for the distance matrices. Example: weights = c(1,2,3). NULL (default): all matrices have same weight in the calculation of W. }
  \item{mult}{ Method for correcting P-values in multiple testing. The methods are "holm" (default), "sidak", and "bonferroni". The Bonferroni correction is overly conservative; it is not recommended. It is included to allow comparisons with the other methods. }
  \item{mantel}{ TRUE: Mantel statistics will be computed from ranked distances, as well as permutational P-values. FALSE (default): Mantel statistics and tests will not be computed. }
  \item{silent}{ TRUE: informative messages will not be printed, but stopping messages will. Option useful for simulation work. FALSE: informative messages will be printed. }
}
\details{
\code{Dmat} must contain two or more distance matrices, listed one after the other, all of the same size, and corresponding to the same objects in the same order. Raw data tables can be transformed into distance matrices before comparison with other such distance matrices, or with data that have been obtained as distance matrices, e.g. serological or DNA hybridization data. The distances will be transformed to ranks before computation of the coefficient of concordance and other statistics.

\code{CADM.global} tests the global null hypothesis that all matrices are incongruent. If the global null is rejected, function \code{CADM.post} can be used to identify the concordant (H0 rejected) and discordant matrices (H0 not rejected) in the group. If a distance matrix has a negative value for the \code{Mantel.mean} statistic, that matrix clearly does not belong to the group. Remove that matrix (if there are more than one, remove first the matrix that has the most strongly negative value for \code{Mantel.mean}) and run the analysis again.

The corrections used for multiple testing are applied to the list of P-values (P) produced in the a posteriori tests; they take into account the number of tests (k) carried out simulatenously (number of matrices, parameter \code{nmat}).

The Holm correction is computed after ordering the P-values in a list with the smallest value to the left. Compute adjusted P-values as:

\deqn{P_{corr} = (k-i+1)*P}{P_corr = (k-i+1)*P}

where i is the position in the ordered list. Final step: from left to right, if an adjusted \eqn{P_{corr}}{P_corr} in the ordered list is smaller than the one occurring at its left, make the smallest one equal to the largest one.

The Sidak correction is:

\deqn{P_{corr} = 1 - (1 - P)^k}{P_corr = 1 - (1 - P)^k}

The Bonferonni correction is:

\deqn{P_{corr} = k*P}{P_corr = k*P}
}

\value{

\code{CADM.global} produces a small table containing the W, Chi2, and Prob.perm statistics described in the following list.
\code{CADM.post} produces a table stored in element \code{A_posteriori_tests}, containing Mantel.mean, Prob, and Corrected.prob statistics in rows; the columns correspond to the k distance matrices under study, labeled Dmat.1 to Dmat.k.
If parameter \code{mantel} is TRUE, tables of Mantel statistics and P-values are computed among the matrices.

  \item{W }{Kendall's coefficient of concordance, W (Kendall and Babington Smith 1939; see also Legendre 2010). }
  \item{Chi2 }{Friedman's chi-square statistic (Friedman 1937) used in the permutation test of W. }
  \item{Prob.perm }{Permutational probability. }

  \item{Mantel.mean }{Mean of the Mantel correlations, computed on rank-transformed distances, between the distance matrix under test and all the other matrices in the study. }
  \item{Prob }{Permutational probabilities, uncorrected. }
  \item{Corrected prob }{Permutational probabilities corrected using the method selected in parameter \code{mult}. }

  \item{Mantel.cor }{Matrix of Mantel correlations, computed on rank-transformed distances, among the distance matrices. }
  \item{Mantel.prob }{One-tailed P-values associated with the Mantel correlations of the previous table. The probabilities are computed in the right-hand tail. H0 is tested against the alternative one-tailed hypothesis that the Mantel correlation under test is positive. No correction is made for multiple testing. }
}

\references{
Campbell, V., Legendre, P. and Lapointe, F.-J. (2009) Assessing congruence among ultrametric distance matrices. \emph{Journal of Classification}, \bold{26}, 103--117.

Campbell, V., Legendre, P. and Lapointe, F.-J. (2011) The performance of the Congruence Among Distance Matrices (CADM) test in phylogenetic analysis. \emph{BMC Evolutionary Biology}, \bold{11}, 64.

Friedman, M. (1937) The use of ranks to avoid the assumption of normality implicit in the analysis of variance. \emph{Journal of the American Statistical Association}, \bold{32}, 675--701.

Kendall, M. G. and Babington Smith, B. (1939) The problem of m rankings. \emph{Annals of Mathematical Statistics}, \bold{10}, 275--287.

Lapointe, F.-J., Kirsch, J. A. W. and Hutcheon, J. M. (1999) Total evidence, consensus, and bat phylogeny: a distance-based approach. \emph{Molecular Phylogenetics and Evolution}, \bold{11}, 55--66.

Legendre, P. (2010) Coefficient of concordance. Pp. 164-169 in: Encyclopedia of Research Design, Vol. 1. N. J. Salkind, ed. SAGE Publications, Inc., Los Angeles.

Legendre, P. and Lapointe, F.-J. (2004) Assessing congruence among distance matrices: single malt Scotch whiskies
revisited. \emph{Australian and New Zealand Journal of Statistics}, \bold{46}, 615--629.

Legendre, P. and Lapointe, F.-J. (2005) Congruence entre matrices de distance. P. 178-181 in: Makarenkov, V., G. Cucumel et F.-J. Lapointe [eds] Comptes rendus des 12emes Rencontres de la Societe Francophone de Classification, Montreal, 30 mai - 1er juin 2005.

Siegel, S. and Castellan, N. J., Jr. (1988) \emph{Nonparametric statistics for the behavioral sciences. 2nd edition}. New York: McGraw-Hill.
}

\author{Pierre Legendre, Universite de Montreal}

\examples{
# Examples 1 and 2: 5 genetic distance matrices computed from simulated DNA
# sequences representing 50 taxa having evolved along additive trees with
# identical evolutionary parameters (GTR+ Gamma + I). Distance matrices were
# computed from the DNA sequence matrices using a p distance corrected with the
# same parameters as those used to simulate the DNA sequences. See Campbell et
# al. (2009) for details.

# Example 1: five independent additive trees. Data provided by V. Campbell.

data(mat5Mrand)
res.global <- CADM.global(mat5Mrand, 5, 50)

# Example 2: three partly similar trees, two independent trees.
# Data provided by V. Campbell.

data(mat5M3ID)
res.global <- CADM.global(mat5M3ID, 5, 50)
res.post   <- CADM.post(mat5M3ID, 5, 50, mantel=TRUE)

# Example 3: three matrices respectively representing Serological
# (asymmetric), DNA hybridization (asymmetric) and Anatomical (symmetric)
# distances among 9 families. Data from Lapointe et al. (1999).

data(mat3)
res.global <- CADM.global(mat3, 3, 9, nperm=999)
res.post   <- CADM.post(mat3, 3, 9, nperm=999, mantel=TRUE)

# Example 4, showing how to bind two D matrices (cophenetic matrices
# in this example) into a file using rbind(), then run the global test.

a <- rtree(5)
b <- rtree(5)
A <- cophenetic(a)
B <- cophenetic(b)
x <- rownames(A)
B <- B[x, x]
M <- rbind(A, B)
CADM.global(M, 2, 5)
}

\keyword{ multivariate }
\keyword{ nonparametric }
