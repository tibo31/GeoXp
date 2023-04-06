\name{histnbmap}
\alias{histnbmap}
\title{Interactive histogram of the distances between two neighbors of a nb object and map}
\description{
The function \code{histnbmap()} draws the histogram of the distances between two neighbors for
a neighbourhood structure given by a \code{nb} object and calculated by \code{nbdists} function
(see \code{spdep} package), and links the graphic with a map.
}

\usage{
histnbmap(sf.obj, nb.obj, longlat = NULL, nbcol = 10,
  type = c("count", "percent", "density"), sup = FALSE, criteria = NULL, 
  carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16, 
  col = "lightblue3", xlab = "", ylab = "count", 
  axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{nb.obj}{object of class nb}
  \item{longlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers;
  if coords is a SpatialPoints object, the value is taken from the object itself}
  \item{nbcol}{number of cells for histogram (10 by default)}
  \item{type}{Character string indicating type of histogram to be drawn. "percent" and "count" give relative frequency and frequency histograms, "density" produces a density scale histogram.}
  \item{sup}{if TRUE, it keeps only the distance of the neighbor the farest }
  \item{criteria}{a vector of boolean of size the number of Spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{"lightblue3" by default, color of bars on the barplot}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
For a selected site j on the map, are represented on the map its neighbours.
For a selected bar on the graph, the corresponding sites are represented on the map with a link which
means that two sites are neighbours.}

\note{
When user select sites on the graph or on the map, he cannot add a selection by using the other graphic.}

\value{
No value returned.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Aragon Y., Thomas-Agnan C., Ruiz-Gazen A., Laurent T.}

\keyword{spatial}

\examples{
##
# data columbus
require("sf")
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
# col.gal.nb is a spatial weight matrix included in spdep package...
library("spdep")
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])

# a simple use of histnbmap
histnbmap(columbus, col.gal.nb, criteria = (columbus$CP == 1),
 xlab = "distance of the neighbor the farest")

\dontrun{ 
require(sf)
if (require(sp, quietly = TRUE)) {
 data(meuse, package = "sp")
 meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
 data(meuse.riv)
 meuse.sr <- st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))),
     crs = 28992)
}

# creation of a spatial weight matrix (class nb) based
# on the Delaunay triangulation
if (require(spdep, quietly = TRUE)) {
  meuse.nb <- tri2nb(st_coordinates(meuse_sf))
}

# a example with some optionswhich shows the limit of
# this kind of spatial weight matrix
histnbmap(meuse_sf, meuse.nb, nbcol = 7,
  carte = meuse.sr)
}
}

\seealso{\code{\link{moranplotmap}}, \code{\link{makeneighborsw}}, \code{\link{normw}}, \code{\link{nonormmoran}} }

