\name{apetools}
\alias{apetools}
\alias{Xplorefiles}
\alias{Xplor}
\alias{editFileExtensions}
\alias{bydir}
\title{Tools to Explore Files}
\description{
  These functions help to find files on the local disk.
}
\usage{
Xplorefiles(from = "HOME", recursive = TRUE, ignore.case = TRUE)
editFileExtensions()
bydir(x)
Xplor(from = "HOME")
}
\arguments{
  \item{from}{the directory where to start the file search; by default,
    the `HOME' directory. Use \code{from = getwd()} to start from the
    current working directory.}
  \item{recursive}{whether to search the subdirectories; \code{TRUE} by
    default.}
  \item{ignore.case}{whether to ignore the case of the file extensions;
    \code{TRUE} by default.}
  \item{x}{a list returned by \code{Xplorefiles}.}
}
\details{
  \code{Xplorefiles} looks for all files with a specified extension in
  their names. The default is to look for the following file types:
  CLUSTAL (.aln), FASTA (.fas, .fasta), FASTQ (.fq, .fastq), NEWICK
  (.nwk, .newick, .tre, .tree), NEXUS (.nex, .nexus), and PHYLIP
  (.phy). This list can be modified with \code{editFileExtensions}.

  \code{bydir} sorts the list of files by directories.

  \code{Xplor} combines the other operations and opens the results in
  a Web browser with clickable links to the directories and files.
}
\value{
  \code{Xplorefiles} returns a list. \code{bydir} prints the file
  listings on the console.
}
\author{Emmanuel Paradis}
\examples{
\dontrun{
x <- Xplorefiles()
x # all data files on your disk
bydir(x) # sorted by directories
bydir(x["fasta"]) # only the FASTA files
Xplorefiles(getwd(), recursive = FALSE) # look only in current dir
Xplor()
}}
\keyword{manip}
