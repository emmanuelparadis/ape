\name{alex}
\alias{alex}
\title{Alignment Explorer With Multiple Devices}
\description{
  This function helps to explore DNA alignments by zooming in. The user
  clicks twice defining the opposite corners of the portion which is
  extracted and drawned on a new window.
}
\usage{
alex(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"}.}
  \item{\dots}{further arguments to pass to \code{image.DNAbin}.}
}
\details{
  This function works with a DNA alignment (freshly) plotted on an
  interactive graphical device (i.e., not a file) with \code{image}.
  After calling \code{alex}, the user clicks twice defining a rectangle
  in the alignment, then this portion of the alignment is extacted and
  plotted on a \emph{new} window. The user can click as many times on
  the alignment. The process is stopped by a right-click. If the user
  clicks twice outside the alignment, a message ``Try again!'' is
  printed.

  Each time \code{alex} is called, the alignment is plotted on a new
  window without closing or deleting those possibly already plotted.

  In all cases, the device where \code{x} is plotted is the active
  window after the operation. It should \emph{not} be closed during the
  whole process.
}
\value{NULL}
\author{Emmanuel Paradis}
\seealso{
  \code{\link{image.DNAbin}}, \code{\link{trex}}, \code{\link{alview}}
}
\examples{
\dontrun{
data(woodmouse)
image(woodmouse)
alex(woodmouse)
}}
\keyword{hplot}
