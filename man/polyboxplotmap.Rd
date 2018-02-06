\name{polyboxplotmap}
\alias{polyboxplotmap}
\title{Interactive polyboxplot and map}
\description{
Parallel Boxplots of a numerical variable by levels of a factor. 
It interacts with a map.
}
\usage{
polyboxplotmap(sp.obj, names.var, varwidth = FALSE, names.arg = "",
 names.attr = names(sp.obj), criteria = NULL, carte = NULL, identify = FALSE, 
 cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "count", 
 axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{names.var}{a vector of 2 characters; first character corrseponds to the name of a factor and second character corresponds to a numeric value}
  \item{varwidth}{if \code{varwidth} is TRUE, the boxes are drawn with widths proportional to the square-roots of the number of observations in the groups}
  \item{names.arg}{a vector of level names (for factor)}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean of size the number of Spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
  \item{cex.lab}{character size of label}
  \item{pch}{a vector of symbol which must be equal to the number of group else all sites are printed in pch[1]}
  \item{col}{a vector of colors which must be equal to the number of group else all sites and all box are printed in col[1]}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
Each site is associated to a value of a level of $names.var[1]$ and of a value of $names.var[2]$.
There is an interactivity only in one direction : the sites selected by quartile-intervals on one of the boxplots 
are then represented on the map in red (or colors according to the options). 
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{dblehistomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{dbledensitymap}}}


\examples{
######
# data eire
require("rgdal")
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])

# example of use of polyboxplotmap
polyboxplotmap(eire, c("pale", "POPCHG"), names.arg = c("Outside Pale", "Pale"), 
 xlab = "Appartenance to the region of Pale", col = c("pink", "violet"))
}
\keyword{spatial}
\keyword{multivariate} 