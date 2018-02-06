\name{boxplotmap}
\alias{boxplotmap}
\title{Interactive boxplot and map}
\description{
The function \code{boxplotmap()} draws a boxplot of the given variable \code{name.var}
and a map with site of coordinates \code{coordinates(sp.obj)}. }

\usage{
boxplotmap(sp.obj, name.var, names.attr = names(sp.obj), criteria = NULL,
 carte = NULL, identify = FALSE, cex.lab = 0.8, pch = 16, col = "lightblue3", 
 xlab = "", ylab = "", axes = FALSE, lablong = "", lablat = "")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean of size the number of spatial sites, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{"lightblue3" by default, color of bars on the boxplot}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
There is an interactivity only in one direction : the sites selected by interquartile on the boxplot 
are represented on the map in red. 
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{histomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}

\examples{
######
# data columbus
require("rgdal")
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])

# a basic use of boxplotmap
boxplotmap(columbus, "CRIME", criteria = (columbus@data$CP == 1), 
 xlab = "Crime", ylab = "Percent", identify = TRUE, cex.lab = 0.7)

\dontrun{ 
# data boston
data(boston)

# creation of a Spatial object
boston.sp <- SpatialPoints(cbind(boston.c$LON, boston.c$LAT))
# ... and then by integrating other variables to create SpatialPointsDataFrame
boston.spdf <- SpatialPointsDataFrame(boston.sp, boston.c)

# a simple use of boxplotmap
boxplotmap(boston.spdf, "MEDV", criteria = (boston.c$CHAS == 1))
}
}

\keyword{spatial}
\keyword{univar} 