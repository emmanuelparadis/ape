\name{parafit}
\alias{parafit}
\alias{print.parafit}
\alias{gopher.D}
\alias{lice.D}
\alias{HP.links}
\title{ Test of host-parasite coevolution }
\description{
Function \code{\link{parafit}} tests the hypothesis of coevolution between a clade of hosts and a clade of parasites. The null hypothesis (H0) of the global test is that the evolution of the two groups, as revealed by the two phylogenetic trees and the set of host-parasite association links, has been independent. Tests of individual host-parasite links are also available as an option.

The method, which is described in detail in Legendre et al. (2002), requires some estimates of the phylogenetic trees or phylogenetic distances, and also a description of the host-parasite associations (H-P links) observed in nature.
}
\usage{
parafit(host.D, para.D, HP, nperm = 999, test.links = FALSE,
        seed = NULL, correction = "none", silent = FALSE)
}

\arguments{
  \item{host.D }{ A matrix of phylogenetic or patristic distances among the hosts (object class: \code{matrix}, \code{data.frame} or \code{dist}). A matrix of patristic distances exactly represents the information in a phylogenetic tree. }
  \item{para.D }{ A matrix of phylogenetic or patristic distances among the parasites (object class: \code{matrix}, \code{data.frame} or \code{dist}). A matrix of patristic distances exactly represents the information in a phylogenetic tree. }
  \item{HP }{ A rectangular matrix with hosts as rows and parasites as columns. The matrix contains 1's when a host-parasite link has been observed in nature between the host in the row and the parasite in the column, and 0's otherwise. }
  \item{nperm}{ Number of permutations for the tests. If \code{nperm =
   0}, permutation tests will not be computed. The default value is \code{nperm = 999}. For large data files, the permutation test is rather slow since the permutation procedure is not compiled. }
  \item{test.links }{ \code{test.links = TRUE} will test the significance of individual host-parasite links. Default: \code{test.links = FALSE}. }
  \item{seed }{ \code{seed = NULL} (default): a seed is chosen at random by the function. That seed is used as the starting point for all tests of significance, i.e. the global H-P test and the tests of individual H-P links if they are requested. Users can select a seed of their choice by giving any integer value to \code{seed}, for example \code{seed = -123456}. Running the function again with the same seed value will produce the exact same test results. }
  \item{correction}{ Correction methods for negative eigenvalues (details below): \code{correction="lingoes"} and \code{correction="cailliez"}. Default value: \code{"none"}.  }
  \item{silent}{ Informative messages and the time to compute the tests will not be written to the \R console if silent=TRUE. Useful when the function is called by a numerical simulation function. }
}

\details{
Two types of test are produced by the program: a global test of coevolution and, optionally, a test on the individual host-parasite (H-P) link.

The function computes principal coordinates for the host and the parasite distance matrices. The principal coordinates (all of them) act as a complete representation of either the phylogenetic distance matrix or the phylogenetic tree.

Phylogenetic distance matrices are normally Euclidean. Patristic distance matrices are additive, thus they are metric and Euclidean. Euclidean matrices are fully represented by real-valued principal coordinate axes. For non-Euclidean matrices, negative eigenvalues are produced; complex principal coordinate axes are associated with the negative eigenvalues. So, the program rejects matrices that are not Euclidean and stops.

Negative eigenvalues can be corrected for by one of two methods: the Lingoes or the Caillez correction. It is up to the user to decide which correction method should be applied. This is done by selecting the option \code{correction="lingoes"} or \code{correction="cailliez"}. Details on these correction methods are given in the help file of the \code{pcoa} function.

The principle of the global test is the following (H0: independent evolution of the hosts and parasites): (1) Compute matrix D = C t(A) B. Note: D is a fourth-corner matrix (sensu Legendre et al. 1997), where A is the H-P link matrix, B is the matrix of principal coordinates computed from the host.D matrix, and C is the matrix of principal coordinates computed from the para.D matrix. (2) Compute the statistic ParaFitGlobal, the sum of squares of all values in matrix D. (3) Permute at random, separately, each row of matrix A, obtaining matrix A.perm. Compute D.perm = C %*% t(A.perm) %*% B, and from it, compute a permuted value ParaFitGlobal.perm for the statistic. Save this value in a vector trace.perm for the tests of individual links (below). (4) Repeat step 4 a large number of times. (5) Add the reference value of ParaFitGlobal to the distribution of ParaFitGlobal.perm values. Calculate the permutational probability associated to ParaFitGlobal.

The test of each individual H-P link is carried out as follows (H0: this particular link is random): (1) Remove one link (k) from matrix A. (2) Compute matrix D = C t(A) B. (3a) Compute trace(k), the sum of squares of all values in matrix D. (3b) Compute the statistic ParaFitLink1 = (trace - trace(k)) where trace is the ParaFitGlobal statistic. (3c) Compute the statistic ParaFitLink2 = (trace - trace(k)) / (tracemax - trace) where tracemax is the maximum value that can be taken by trace. (4) Permute at random, separately, each row of matrix A, obtaining A.perm. Use the same sequences of permutations as were used in the test of ParaFitGlobal. Using the values of trace and trace.perm saved during the global test, compute the permuted values of the two statistics, ParaFit1.perm and ParaFit2.perm. (5) Repeat step 4 a large number of times. (6) Add the reference value of ParaFit1 to the distribution of ParaFit1.perm values; add the reference value of ParaFit2 to the distribution of ParaFit2.perm values. Calculate the permutational probabilities associated to ParaFit1 and ParaFit2.

The \code{print.parafit} function prints out the results of the global test and, optionally, the results of the tests of the individual host-parasite links.
}

\value{

  \item{ParaFitGlobal }{The statistic of the global H-P test. }
  \item{p.global }{The permutational p-value associated with the ParaFitGlobal statistic. }
  \item{link.table }{The results of the tests of individual H-P links, including the ParaFitLink1 and ParaFitLink2 statistics and the p-values obtained from their respective permutational tests. }
  \item{para.per.host }{Number of parasites per host. }
  \item{host.per.para }{Number of hosts per parasite. }
  \item{nperm }{Number of permutations for the tests. }
  }

\author{ Pierre Legendre, Universite de Montreal }

\references{
Hafner, M. S, P. D. Sudman, F. X. Villablanca, T. A. Spradling, J. W. Demastes and S. A. Nadler. 1994. Disparate rates of molecular evolution in cospeciating hosts and parasites. \emph{Science}, \bold{265}, 1087--1090.

Legendre, P., Y. Desdevises and E. Bazin. 2002. A statistical test for host-parasite coevolution. \emph{Systematic Biology}, \bold{51(2)}, 217--234.
}

\seealso{\code{\link{pcoa}} }

\examples{
## Gopher and lice data from Hafner et al. (1994)

data(gopher.D)
data(lice.D)
data(HP.links)

res <- parafit(gopher.D, lice.D, HP.links, nperm=99, test.links=TRUE)
# res     # or else: print(res)
}

\keyword{ multivariate }
