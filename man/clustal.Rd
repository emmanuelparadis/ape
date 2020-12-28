\name{clustal}
\alias{clustal}
\alias{clustalomega}
\alias{muscle}
\alias{tcoffee}
\title{Multiple Sequence Alignment with External Applications}
\description{
  These functions call their respective program from \R to align a set
  of nucleotide sequences of class \code{"DNAbin"} or
  \code{"AAbin"}. The application(s) must be installed seperately and it
  is highly recommended to do this so that the executables are in a
  directory located on the PATH of the system.
}
\usage{
clustal(x, y, guide.tree, pw.gapopen = 10, pw.gapext = 0.1,
        gapopen = 10, gapext = 0.2, exec = NULL, MoreArgs = "",
        quiet = TRUE, original.ordering = TRUE, file)
clustalomega(x, y, guide.tree, exec = NULL,MoreArgs = "",
              quiet = TRUE, original.ordering = TRUE, file)
muscle(x, y, guide.tree, exec, MoreArgs = "",
        quiet = TRUE, original.ordering = TRUE, file)
tcoffee(x, exec = "t_coffee", MoreArgs = "", quiet = TRUE,
        original.ordering = TRUE)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"} or \code{"AAbin"} (can be
    missing).}
  \item{y}{an object of class \code{"DNAbin"} or \code{"AAbin"} used for
    profile alignment (can be missing).}
  \item{guide.tree}{guide tree, an object of class \code{"phylo"} (can
    be missing).}
  \item{pw.gapopen, pw.gapext}{gap opening and gap extension penalties
    used by Clustal during pairwise alignments.}
  \item{gapopen, gapext}{idem for global alignment.}
  \item{exec}{a character string giving the name of the program, with
    its path if necessary. \code{clustal} tries to guess this argument
    depending on the operating system (see details).}
  \item{MoreArgs}{a character string giving additional options.}
  \item{quiet}{a logical: the default is to not print on \R's console the
    messages from the external program.}
  \item{original.ordering}{a logical specifying whether to return the
    aligned sequences in the same order than in \code{x} (\code{TRUE} by
    default).}
  \item{file}{a file with its path if results should be stored (can be
    missing).}
}
\details{
  It is highly recommended to install the executables properly so that
  they are in a directory located on the PATH (i.e., accessible from any
  other directory). Alternatively, the full path to the executable
  may be given (e.g., \code{exec = "~/muscle/muscle"}), or a (symbolic)
  link may be copied in the working directory. For Debian and its
  derivatives (e.g., Ubuntu), it is recommended to use the binaries
  distributed by Debian.

  \code{clustal} tries to guess the name of the executable program
  depending on the operating system. Specifically, the followings are
  used: ``clustalw'' under Linux, ``clustalw2'' under MacOS, and
  ``clustalw2.exe'' under Windows. For \code{clustalomega},
  ``clustalo[.exe]'' is the default on all systems (with no specific
  path).

  When called without arguments (i.e., \code{clustal()}, \dots), the
  function prints the options of the program which may be passed to
  \code{MoreArgs}.

  Since \pkg{ape} 5.1, \code{clustal}, \code{clustalomega}, and
  \code{muscle} can align AA sequences as well as DNA sequences.
}
\value{
  an object of class \code{"DNAbin"} or \code{"AAbin"} with the aligned
  sequences.
}
\references{
  Chenna, R., Sugawara, H., Koike, T., Lopez, R., Gibson, T. J.,
  Higgins, D. G. and Thompson, J. D. (2003) Multiple sequence alignment
  with the Clustal series of programs. \emph{Nucleic Acids Research}
  \bold{31}, 3497--3500. \url{http://www.clustal.org/}

  Edgar, R. C. (2004) MUSCLE: Multiple sequence alignment with high
  accuracy and high throughput. \emph{Nucleic Acids Research},
  \bold{32}, 1792--1797.
  \url{http://www.drive5.com/muscle/muscle_userguide3.8.html}

  Notredame, C., Higgins, D. and Heringa, J. (2000) T-Coffee: A novel
  method for multiple sequence alignments. \emph{Journal of Molecular
  Biology}, \bold{302}, 205--217.
  \url{http://www.tcoffee.org/}

  Sievers, F., Wilm, A., Dineen, D., Gibson, T. J., Karplus, K., Li, W.,
  Lopez, R., McWilliam, H., Remmert, M., S\"oding, J., Thompson,
  J. D. and Higgins, D. G. (2011) Fast, scalable generation of
  high-quality protein multiple sequence alignments using Clustal
  Omega. \emph{Molecular Systems Biology}, \bold{7}, 539.
  \url{http://www.clustal.org/}
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{image.DNAbin}}, \code{\link{del.gaps}},
  \code{\link{all.equal.DNAbin}}, \code{\link{alex}},
  \code{\link{alview}}, \code{\link{checkAlignment}}
}
\examples{
\dontrun{
### display the options:
clustal()
clustalomega()
muscle()
tcoffee()

data(woodmouse)
### open gaps more easily:
clustal(woodmouse, pw.gapopen = 1, pw.gapext = 1)
### T-Coffee requires negative values (quite slow; muscle() is much faster):
tcoffee(woodmouse,  MoreArgs = "-gapopen=-10 -gapext=-2")
}}
\keyword{manip}
