\name{dist.dna}
\alias{dist.dna}
\title{Pairwise Distances from DNA Sequences}
\usage{
dist.dna(x, model = "K80", variance = FALSE,
         gamma = FALSE, pairwise.deletion = FALSE,
         base.freq = NULL, as.matrix = FALSE)
}
\arguments{
  \item{x}{a matrix or a list containing the DNA sequences; this must be
    of class \code{"DNAbin"} (use \code{\link{as.DNAbin}} is they are
    stored as character).}
  \item{model}{a character string specifying the evolutionary model to be
    used; must be one of \code{"raw"}, \code{"N"}, \code{"TS"},
    \code{"TV"}, \code{"JC69"}, \code{"K80"} (the default),
    \code{"F81"}, \code{"K81"}, \code{"F84"}, \code{"BH87"},
    \code{"T92"}, \code{"TN93"}, \code{"GG95"}, \code{"logdet"},
    \code{"paralin"}, \code{"indel"}, or \code{"indelblock"}.}
  \item{variance}{a logical indicating whether to compute the variances
    of the distances; defaults to \code{FALSE} so the variances are not
    computed.}
  \item{gamma}{a value for the gamma parameter possibly used to apply a
    correction to the distances (by default no correction is applied).}
  \item{pairwise.deletion}{a logical indicating whether to delete the
    sites with missing data in a pairwise way. The default is to delete
    the sites with at least one missing data for all sequences (ignored
    if \code{model = "indel"} or \code{"indelblock"}).}
  \item{base.freq}{the base frequencies to be used in the computations
    (if applicable). By default, the base frequencies are computed from
    the whole set of sequences.}
  \item{as.matrix}{a logical indicating whether to return the results as
    a matrix. The default is to return an object of class
    \link[stats]{dist}.}
}
\description{
  This function computes a matrix of pairwise distances from DNA
  sequences using a model of DNA evolution. Eleven substitution models
  (and the raw distance) are currently available.
}
\details{
  The molecular evolutionary models available through the option
  \code{model} have been extensively described in the literature. A
  brief description is given below; more details can be found in the
  references.

\itemize{
  \item{\code{raw}, \code{N}: }{This is simply the proportion or the number of
    sites that differ between each pair of sequences. This may be useful
    to draw ``saturation plots''. The options \code{variance} and
    \code{gamma} have no effect, but \code{pairwise.deletion} can.}

  \item{\code{TS}, \code{TV}: }{These are the numbers of transitions and
    transversions, respectively.}

  \item{\code{JC69}: }{This model was developed by Jukes and Cantor (1969). It
    assumes that all substitutions (i.e. a change of a base by another
    one) have the same probability. This probability is the same for all
    sites along the DNA sequence. This last assumption can be relaxed by
    assuming that the substition rate varies among site following a
    gamma distribution which parameter must be given by the user. By
    default, no gamma correction is applied. Another assumption is that
    the base frequencies are balanced and thus equal to 0.25.}

  \item{\code{K80}: }{The distance derived by Kimura (1980), sometimes referred
    to as ``Kimura's 2-parameters distance'', has the same underlying
    assumptions than the Jukes--Cantor distance except that two kinds of
    substitutions are considered: transitions (A <-> G, C <-> T), and
    transversions (A <-> C, A <-> T, C <-> G, G <-> T). They are assumed
    to have different probabilities. A transition is the substitution of
    a purine (C, T) by another one, or the substitution of a pyrimidine
    (A, G) by another one. A transversion is the substitution of a
    purine by a pyrimidine, or vice-versa. Both transition and
    transversion rates are the same for all sites along the DNA
    sequence. Jin and Nei (1990) modified the Kimura model to allow for
    variation among sites following a gamma distribution. Like for the
    Jukes--Cantor model, the gamma parameter must be given by the
    user. By default, no gamma correction is applied.}

  \item{\code{F81}: }{Felsenstein (1981) generalized the Jukes--Cantor model
    by relaxing the assumption of equal base frequencies. The formulae
    used in this function were taken from McGuire et al. (1999)}.

  \item{\code{K81}: }{Kimura (1981) generalized his model (Kimura 1980) by
    assuming different rates for two kinds of transversions: A <-> C and
    G <-> T on one side, and A <-> T and C <-> G on the other. This is
    what Kimura called his ``three substitution types model'' (3ST), and
    is sometimes referred to as ``Kimura's 3-parameters distance''}.

  \item{\code{F84}: }{This model generalizes K80 by relaxing the assumption
    of equal base frequencies. It was first introduced by Felsenstein in
    1984 in Phylip, and is fully described by Felsenstein and Churchill
    (1996). The formulae used in this function were taken from McGuire
    et al. (1999)}.

  \item{\code{BH87}: }{Barry and Hartigan (1987) developed a distance based
    on the observed proportions of changes among the four bases. This
    distance is not symmetric.}

  \item{\code{T92}: }{Tamura (1992) generalized the Kimura model by relaxing
    the assumption of equal base frequencies. This is done by taking
    into account the bias in G+C content in the sequences. The
    substitution rates are assumed to be the same for all sites along
    the DNA sequence.}

  \item{\code{TN93}: }{Tamura and Nei (1993) developed a model which assumes
    distinct rates for both kinds of transition (A <-> G versus C <->
    T), and transversions. The base frequencies are not assumed to be
    equal and are estimated from the data. A gamma correction of the
    inter-site variation in substitution rates is possible.}

  \item{\code{GG95}: }{Galtier and Gouy (1995) introduced a model where the
    G+C content may change through time. Different rates are assumed for
    transitons and transversions.}

  \item{\code{logdet}: }{The Log-Det distance, developed by Lockhart et
    al. (1994), is related to BH87. However, this distance is
    symmetric. Formulae from Gu and Li (1996) are used.
    \code{dist.logdet} in \pkg{phangorn} uses a different
    implementation that gives substantially different distances for
    low-diverging sequences.}

  \item{\code{paralin}: }{Lake (1994) developed the paralinear distance which
    can be viewed as another variant of the Barry--Hartigan distance.}

  \item{\code{indel}: }{this counts the number of sites where there is an
    insertion/deletion gap in one sequence and not in the other.}

  \item{\code{indelblock}: }{same than before but contiguous gaps are
    counted as a single unit. Note that the distance between \code{-A-} and
    \code{A--} is 3 because there are three different blocks of gaps, whereas
    the ``indel'' distance will be 2.}
}}
\note{
  If the sequences are very different, most evolutionary distances are
  undefined and a non-finite value (Inf or NaN) is returned. You may do
  \code{dist.dna(, model = "raw")} to check whether some values are
  higher than 0.75.
}
\value{
  an object of class \link[stats]{dist} (by default), or a numeric
  matrix if \code{as.matrix = TRUE}. If \code{model = "BH87"}, a numeric
  matrix is returned because the Barry--Hartigan distance is not
  symmetric.

  If \code{variance = TRUE} an attribute called \code{"variance"} is
  given to the returned object.
}
\references{
  Barry, D. and Hartigan, J. A. (1987) Asynchronous distance between
  homologous DNA sequences. \emph{Biometrics}, \bold{43}, 261--276.

  Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a
  maximum likelihood approach. \emph{Journal of Molecular Evolution},
  \bold{17}, 368--376.

  Felsenstein, J. and Churchill, G. A. (1996) A Hidden Markov model
  approach to variation among sites in rate of evolution.
  \emph{Molecular Biology and Evolution}, \bold{13}, 93--104.

  Galtier, N. and Gouy, M. (1995) Inferring phylogenies from DNA
  sequences of unequal base compositions. \emph{Proceedings of the
    National Academy of Sciences USA}, \bold{92}, 11317--11321.

  Gu, X. and Li, W.-H. (1996) Bias-corrected paralinear and LogDet
  distances and tests of molecular clocks and phylogenies under
  nonstationary nucleotide frequencies. \emph{Molecular Biology and
    Evolution}, \bold{13}, 1375--1383.

  Jukes, T. H. and Cantor, C. R. (1969) Evolution of protein
  molecules. in \emph{Mammalian Protein Metabolism}, ed. Munro, H. N.,
  pp. 21--132, New York: Academic Press.

  Kimura, M. (1980) A simple method for estimating evolutionary rates of
  base substitutions through comparative studies of nucleotide
  sequences. \emph{Journal of Molecular Evolution}, \bold{16}, 111--120.

  Kimura, M. (1981) Estimation of evolutionary distances between
  homologous nucleotide sequences. \emph{Proceedings of the National
    Academy of Sciences USA}, \bold{78}, 454--458.

  Jin, L. and Nei, M. (1990) Limitations of the evolutionary parsimony
  method of phylogenetic analysis. \emph{Molecular Biology and
    Evolution}, \bold{7}, 82--102.

  Lake, J. A. (1994) Reconstructing evolutionary trees from DNA and
  protein sequences: paralinear distances. \emph{Proceedings of the
    National Academy of Sciences USA}, \bold{91}, 1455--1459.

  Lockhart, P. J., Steel, M. A., Hendy, M. D. and Penny, D. (1994)
  Recovering evolutionary trees under a more realistic model of sequence
  evolution. \emph{Molecular Biology and Evolution}, \bold{11},
  605--602.

  McGuire, G., Prentice, M. J. and Wright, F. (1999). Improved error
  bounds for genetic distances from DNA sequences. \emph{Biometrics},
  \bold{55}, 1064--1070.

  Tamura, K. (1992) Estimation of the number of nucleotide substitutions
  when there are strong transition-transversion and G + C-content
  biases. \emph{Molecular Biology and Evolution}, \bold{9}, 678--687.

  Tamura, K. and Nei, M. (1993) Estimation of the number of nucleotide
  substitutions in the control region of mitochondrial DNA in humans and
  chimpanzees. \emph{Molecular Biology and Evolution}, \bold{10}, 512--526.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{read.GenBank}}, \code{\link{read.dna}},
  \code{\link{write.dna}},  \code{\link{DNAbin}},
  \code{\link{dist.gene}}, \code{\link{cophenetic.phylo}},
  \code{\link[stats]{dist}}
}
\keyword{manip}
\keyword{multivariate}
\keyword{cluster}
