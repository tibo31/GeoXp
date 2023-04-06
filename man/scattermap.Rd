\name{scattermap}
\alias{scattermap}

\title{Interactive scatterplot and map}

\description{
The function \code{scattermap} draws a scatterplot of the given variables indicated
in \code{names.var} and a map with sites
of coordinates \code{coordinates(sf.obj)}. Boxplots of each variable
\code{names.var[1]} and \code{names.var[2]} are represented below the x-axis and y-axis.}

\usage{
scattermap(sf.obj, names.var, lin.reg = TRUE, quantiles = TRUE,
  criteria = NULL, carte = NULL, identify = NULL, 
  cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", 
  axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{names.var}{a vector of two characters : 1st name corresponds to the x-variable,
  2nd name corresponds to the y-variable}
  \item{lin.reg}{If TRUE, drawing of the linear predictor for `response' in linear model}
  \item{quantiles}{a boolean to represent the Additive Quantile Regression Smoothing}
  \item{criteria}{a vector of boolean of size the number of Spatial Units, which permit to represent preselected sites with a cross, using the tcltk window}
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
\details{Each site is associated to a value of \code{names.var[1]} and a value of
\code{names.var[2]}. There is an interactivity between the two windows: the sites selected
by `point' or `polygon' on the scatterplot are represented on the map in red; sites selected
on the map are then represented in red on the scatterplot. Users have
the possibility to draw linear predictor for `response' in linear model (option \code{lin.reg})
 or conditionnal quantile regression spline (option \code{quantiles}).}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{dblehistomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}

\examples{
######
# data on price indices of real estate in France
# data on price indices of real estate in France
data(immob, package = "GeoXp")

# immob is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
require(sf)
immob.sf <- st_as_sf(immob, coords = c("longitude", "latitude"))
# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
midiP <- st_read(system.file("shapes/region.shp", package="GeoXp")[1])

# a example of use
scattermap(immob.sf, c("prix.vente", "prix.location"),
  carte = midiP, xlab = "Average sell price", ylab = "Average rent price",
  identify = "Nom", cex.lab = 0.6)

\dontrun{ 
# data eire
eire <- sf::st_read(system.file("shapes/eire.shp", package="spData")[1])

# example of use of scattermap
scattermap(eire, c("ROADACC", "OWNCONS"), lin.reg = TRUE, 
  xlab = "Roads", ylab = "Consomation Rate", col = "purple")
}
}

\keyword{regression}
\keyword{spatial}
\keyword{multivariate} 