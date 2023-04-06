\name{plot3dmap}
\alias{plot3dmap}

\title{Interactive Plot3d and map}
\description{
The function \code{plot3dmap()} draws a 3d-plot of three given variables $names.var$
and a map with sites of coordinates \code{coordinates(sp.obj)}.
}
\usage{
plot3dmap(sf.obj, names.var, box = TRUE,  
  criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8,
  pch = 16, col = "lightblue3", xlab = "", ylab = "", zlab = "", 
  axes = FALSE, lablong = "", lablat = "")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{names.var}{a vector of three characters; attribute names or column numbers in attribute table}
  \item{box}{a boolean with TRUE for drawing a box on the scatterplot 3d}
  \item{criteria}{a vector of boolean of size the number of Spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{"lightblue3" by default, color of bars on the histogram}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{zlab}{a title for the graphic z-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
Sites selected on the map by `points' or `polygon' are
represented in red in the 3-d plot.
}
\note{
This function uses the \code{rgl} package and open a \code{rgl} device.
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y.,  Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{scattermap}}}
\examples{
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
# an example of plot3dmap
plot3dmap(immob.sf, c("prix.vente", "prix.location", "variation.vente"),
  box = FALSE, carte = midiP, identify = "Nom", cex.lab = 0.5,
  xlab = "prix.vente", ylab = "prix.location", zlab = "variation.vente")

######
# data eire
eire <- st_read(system.file("shapes/eire.shp", package="spData")[1])

# an example of use
plot3dmap(eire, c("A", "RETSALE", "INCOME"), 
  xlab = "A", ylab = "RETSALE", zlab = "INCOME")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
\keyword{multivariate}
