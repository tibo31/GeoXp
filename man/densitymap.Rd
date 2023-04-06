\name{densitymap}
\alias{densitymap}

\title{Kernel density estimates and map}
\description{
The function \code{densitymap} draws kernel density estimates of the variable \code{name.var} with
\code{bkde} and a map with sites of coordinate \code{sf.obj}. Each site is associated
to a value of \code{name.var} and there is interactivity between the two windows.
}

\usage{
densitymap(sf.obj, name.var, kernel = "triweight", 
  criteria = NULL, carte = NULL, identify = "",
  cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", axes = FALSE,
  lablong = "", lablat = "")
}

\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{kernel}{Smoothing kernel (see help(bkde) for list of options)}
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
The user can choose an interval on the density curve by mouse clicking on the graph on the extremities of interval or by specifying directly values. The sites selected by an interval are then represented on the map in red.
The selection of sites on the map by `points' or `polygon' results in the drawing of the kernel densities of the subdistributions corresponding to this subset of sites.
Finally, the user can modify the bandwith parameter with a cursor in the Tk
window (parameter \eqn{\alpha}{alpha}). \eqn{\alpha}{alpha} is the smoothing parameter for the kernel smooth : it represents the mean percentage of sample points involved in the local averaging (example : \eqn{\alpha=20}{alpha=20} means that on average, \eqn{n \times 0.2}{n x 0.2} points are in any interval of length 2h where h is the usual bandwidth).}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. New York: Springer. \cr \cr

Wand M.P. et Jones M.C. (1995), \emph{Kernel Smoothing}, Chapman \& Hall.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{histomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}

\examples{

######
# data columbus
require("sf")
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1])

# columbus is included in the Spatial-Class object
# a very simple use of histomap :
densitymap(columbus, "CRIME")

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


# A basic call of densitymap function
densitymap(immob.sf, "prix.vente", carte = midiP, identify = "Nom",
   xlab = "housing price by square meter", cex.lab = 0.6)
}
}

\keyword{spatial}
\keyword{smooth}
\keyword{univar}
\keyword{multivariate}
