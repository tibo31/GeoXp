\name{neighbourmap}
\alias{neighbourmap}
\title{Neighbour plot and map}
\description{
The function \code{neighbourmap()} identifies spatial outliers by comparing a variable value for
a particular site with these of its neighbouring sites. It draws a scatterplot of the values of the
variable at neighbouring sites for a neighbourhood structure given by a binary weight matrix W and
links this scatterplot with a map.
}

\usage{
neighbourmap(sf.obj, name.var, nb.obj, lin.reg = TRUE, 
  criteria = NULL, carte = NULL, identify = NULL, 
  cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", 
  axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{nb.obj}{object of class nb}
  \item{lin.reg}{If TRUE, drawing the linear curve y=x}
  \item{criteria}{a vector of boolean of size the number of spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
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
For a selected site j on the map, are represented on the map its neighbors, and on the graph, on the x-axis the value of \code{var} for this site, and in y-axis, the values of \code{var} for the neighbouring sites of j. For a selected point on the graph, the corresponding pair of sites is represented on the map with a link.}

\note{
When user selects sites on the graph or on the map, he cannot add a selection by using the other graphic.}

\value{
In the case where user click on \code{save results} button,
a matrix of integer is created as a global variable in \code{last.select} object.
It corresponds to the numbers of spatial unit corresponding to couple of sites selected
just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Aragon Y., Thomas-Agnan C., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}

\examples{
###
# columbus
require("sf")
columbus <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])

# example of use of neighbourmap
neighbourmap(columbus, "CRIME", col.gal.nb,
  criteria = (columbus$CRIME > mean(columbus$CRIME)))

###
# data immob
# data on price indices of real estate in France
data(immob, package = "GeoXp")

# immob is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
require(sf)
immob.sf <- st_as_sf(immob, coords = c("longitude", "latitude"))
# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
midiP <- st_read(system.file("shapes/region.shp", package="GeoXp")[1])

# A spatial weight matrix based on triangulation Delaunay
W.nb <- spdep::tri2nb(cbind(immob$longitude, immob$latitude))

# example of use of neighbourmap
neighbourmap(immob.sf, "prix.vente", W.nb, identify = "Nom", cex.lab = 0.5,
   carte = midiP)
}

\seealso{\code{\link{moranplotmap}},\code{\link{makeneighborsw}},\code{\link{normw}},\code{\link{nonormmoran}} }

