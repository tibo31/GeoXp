\name{dblehistomap}
\alias{dblehistomap}
\title{Double Interactive Histogram and map}
\description{
The function \code{dblehistomap} draws two histograms of the given variables \code{names.var[1]} and
\code{names.var[2]} and a map with sites of coordinates \code{sf.obj}. Each site is associated
to a value of \code{names.var[1]} and \code{names.var[2]} and there is interactivity between the two
windows created.}

\usage{
dblehistomap(sf.obj, names.var, nbcol = c(10, 10), 
  type = c("count", "percent", "density"), 
  criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16,
  col = c("grey", "lightblue3"), xlab = c("", ""), ylab = c("count", "count"),
  axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{names.var}{a vector of 2 characters; attribute name or column number in attribute table}
  \item{nbcol}{a vector of integer indicating number of cells for histogram 1 and histogram 2 (10 for each by default)}
  \item{type}{Character string indicating type of histogram to be drawn. "percent" and "count" give relative frequency and frequency histograms, "density" produces a density scale histogram.}
  \item{criteria}{a vector of boolean of size the number of spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{a vector of character indicating colors for the bars on the histogram 1 and histogram 2}
  \item{xlab}{a vector of character containing titles for the two graphics x-axis}
  \item{ylab}{a vector of character containing titles for the two graphics y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{The selection of sites by `points' or `polygons' on the map results in the
drawing of the red histograms of the subdistributions corresponding to this subset of sites.}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{dblehistomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}

\examples{
######
# data on price indices of real estate in France
data(immob, package = "GeoXp")

# immob is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
require(sf)
immob.sf <- st_as_sf(immob, coords = c("longitude", "latitude"))
# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
midiP <- st_read(system.file("shapes/region.shp", package="GeoXp")[1])

# A basic call of dblehistomap function
dblehistomap(immob.sf, c("prix.vente", "prix.location"),
  carte = midiP, identify = "Nom", cex.lab = 0.6)

\dontrun{ 
# data colombus
require("sf")
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])

# example of use with many options
dblehistomap(columbus, c("HOVAL", "CRIME"), nbcol = c(5, 10), type = "percent", 
  xlab = c("hoval", "crime"), identify = "POLYID", cex.lab = 0.7, ylab = c("percent", "percent"),
  col = c("pink", "orange"), pch = 14, axes = TRUE)
}
}
\keyword{spatial}
\keyword{multivariate}
