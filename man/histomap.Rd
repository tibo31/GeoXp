\name{histomap}
\alias{histomap}

\title{Interactive Histogram and map}

\description{
The function \code{histomap()} draws a histogram of a given variable \code{name.var}
and a map with sites of coordinates \code{coordinates(sf.obj)}. Each site is associated to a value
of \code{name.var} and there is interactivity between the two windows.
}

\usage{
histomap(sf.obj, name.var, nbcol = 10, type = c("count", "percent", "density"),
   criteria = NULL, carte = NULL, identify = NULL, 
   cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", axes = FALSE, 
   lablong = "", lablat = "")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{nbcol}{number of cells for histogram (10 by default)}
  \item{type}{Character string indicating type of histogram to be drawn. "percent" and "count" give relative frequency and frequency histograms, "density" produces a density scale histogram.}
  \item{criteria}{a vector of size n of boolean which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{"lightblue3" by default, color of bars on the histogram}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
Sites selected by a bar on the histogram are represented on the map in red and the values
of sites selected on the map by `points' or `polygon' are
represented in red as a sub-histogram on the histogram.
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}
\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan C., Aragon Y.,  Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{histomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}
\examples{
######
# data columbus
require("sf")
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])

# columbus is included in the Spatial-Class object
# a very simple use of histomap :
histomap(columbus, "CRIME")

\dontrun{ 
# data on price indices of real estate in France
data(immob, package = "GeoXp")

# immob is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
require(sf)
immob.sf <- st_as_sf(immob, coords = c("longitude", "latitude"))
# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
midiP <- st_read(system.file("shapes/region.shp", package="GeoXp")[1])

# A basic call of histomap function
histomap(immob.sf, "prix.vente", carte = midiP, 
  identify = "Nom", cex.lab=0.6)

# ... with all options
histomap(immob.sf, "prix.vente", nbcol = 15, type = "percent",
  criteria = immob$rentabilite > 5, carte = midiP, identify = "Nom", 
  cex.lab = 0.5, pch = 12, 
  col = "pink", xlab = "variation price", ylab = "percent", 
  axes = TRUE, lablong = "x", lablat = "y")
}
}
\keyword{spatial}
\keyword{univar}
