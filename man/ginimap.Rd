\name{ginimap}
\alias{ginimap}
\title{Lorentz curve and map}
\description{
The function \code{ginimap()} computes a Lorentz curve from \code{name.var}
and calculates the Gini Index associated to \code{name.var}.
}
\usage{
ginimap(sp.obj, name.var, names.attr = names(sp.obj), criteria = NULL,
 carte = NULL, identify = FALSE, cex.lab = 0.8, pch = 16, col = "lightblue3",
 xlab = "", ylab = "", axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean of size the number of Spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
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
Users have the possibility to choose a threshold by mouse clicking on the Lorentz curve or by specifying it in the menu. The corresponding
pair (F,G) and the value of \code{name.var} are then printed on the graph and the sites
with a value of \code{name.var} lower or equal to the threshold are then selected on the map.
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\note{The Gini Index is given in the tcltk window (see function \code{gini} for the formula used to calculate it).}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y.,  Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{gini}}}

\examples{
######
# data eire
require("rgdal")
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])
# a basic usage ...
ginimap(eire, "INCOME")

\dontrun{ 
# ... with options
ginimap(eire, "INCOME", criteria = (eire$pale == 1),
 identify = TRUE, pch = 5, col = "orange", 
 axes = TRUE, lablong = "X", lablat = "Y")
}
}


\keyword{spatial}
\keyword{univar}