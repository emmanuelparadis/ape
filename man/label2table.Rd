\name{label2table}
\alias{label2table}
\alias{stripLabel}
\alias{abbreviateGenus}
\title{Label Management}
\description{
  These functions work on a vector of character strings storing bi- or trinomial species names, typically ``Genus_species_subspecies''.
}
\usage{
label2table(x, sep = NULL, as.is = FALSE)
stripLabel(x, species = FALSE, subsp = TRUE, sep = NULL)
abbreviateGenus(x, genus = TRUE, species = FALSE, sep = NULL)
}
\arguments{
  \item{x}{a vector of mode character.}
  \item{sep}{the separator (a single character) between the taxonomic levels (see details).}
  \item{as.is}{a logical specifying whether to convert characters into factors (like in \code{\link[utils]{read.table}}).}
  \item{species, subsp, genus}{a logical specifying whether the taxonomic level is concerned by the operation.}
}
\details{
\code{label2table} returns a data frame with three columns named ``genus'', ``species'', and ``subspecies'' (with \code{NA} if the level is missing).

\code{stripLabel} deletes the subspecies names from the input. If \code{species = TRUE}, the species names are also removed, thus returning only the genus names.

\code{abbreviateGenus} abbreviates the genus names keeping only the first letter. If \code{species = TRUE}, the species names are abbreviated.

By default, these functions try to guess what is the separator between the genus, species and/or subspecies names. If an underscore is present in the input, then this character is assumed to be the separator; otherwise, a space. If this does not work, you can specify \code{sep} to its appropriate value.
}
\value{
  A vector of mode character or a data frame.
}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{makeLabel}}, \code{\link{makeNodeLabel}},
  \code{\link{mixedFontLabel}}, \code{\link{updateLabel}},
  \code{\link{checkLabel}}
}
\examples{
x <- c("Panthera_leo", "Panthera_pardus", "Panthera_onca", "Panthera_uncia",
       "Panthera_tigris_altaica", "Panthera_tigris_amoyensis")
label2table(x)
stripLabel(x)
stripLabel(x, TRUE)
abbreviateGenus(x)
abbreviateGenus(x, species = TRUE)
abbreviateGenus(x, genus = FALSE, species = TRUE)
}
\keyword{manip}
