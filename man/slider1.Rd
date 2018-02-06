\name{slider1}
\alias{slider1}
\title{Scroll bar}

\description{
Creates a scroll bar for modifying the value of a parameter.
}
\usage{
slider1(fenetre, refresh.code, names.slide, minima, maxima, resolutions,
starts, no=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fenetre}{number of windows}
  \item{refresh.code}{name of function called after modifying the value of a parameter}
  \item{names.slide}{title for scroll bar}
  \item{minima}{minimum value of parameter}
  \item{maxima}{maximum value of parameter}
  \item{resolutions}{scale}
  \item{starts}{Initial Value}
  \item{no}{number of scroll bar}
}
\details{
This function is used in the functions which draw a curve to modify the smoothing parameter. 
This is not an interactive function}

\value{Draws a Tk window.}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}
\seealso{\code{\link{densitymap}},\code{\link{dbledensitymap}}}
