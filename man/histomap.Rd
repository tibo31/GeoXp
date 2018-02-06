\name{histomap}
\alias{histomap}

\title{Interactive Histogram and map}

\description{
The function \code{histomap()} draws a histogram of a given variable \code{name.var}
and a map with sites of coordinates \code{coordinates(sp.obj)}. Each site is associated to a value
of \code{name.var} and there is interactivity between the two windows.
}

\usage{
histomap(sp.obj, name.var, nbcol = 10, type = c("count", "percent", "density"),
 names.attr = names(sp.obj), criteria = NULL, carte = NULL, identify = FALSE, 
 cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", axes = FALSE, 
 lablong = "", lablat = "")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{nbcol}{number of cells for histogram (10 by default)}
  \item{type}{Character string indicating type of histogram to be drawn. "percent" and "count" give relative frequency and frequency histograms, "density" produces a density scale histogram.}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of size n of boolean which permit to represent preselected sites with a cross, using the tcltk window}
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
# data oldcol
require("rgdal")
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])

# columbus is included in the Spatial-Class object
# a very simple use of histomap :
histomap(columbus, "CRIME")

\dontrun{ 
# data on price indices of real estate in France
data(immob)

# immob is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
# to create Spatial Points object ...
immob.sp <- SpatialPoints(cbind(immob$longitude,immob$latitude))
# ... and then by integrating other variables to create SpatialPointsDataFrame
immob.spdf <- SpatialPointsDataFrame(immob.sp, immob)
# we just give names to the spatial units...
row.names(immob.spdf) <- immob$Nom
# For more details, see vignette('sp', package="sp")

# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
midiP <- readOGR(system.file("shapes/region.shp", package="GeoXp")[1])
cont_midiP <- spdf2list(midiP[-c(22,23),])$poly

# A basic call of histomap function
histomap(immob.spdf, "prix.vente", carte = cont_midiP, 
 identify=TRUE, cex.lab=0.6)

# ... with all options
histomap(immob.spdf,7, nbcol = 15, type = "percent",
 names.attr = names(immob), criteria = immob$rentabilite > 5, 
 carte = cont_midiP, identify = TRUE, cex.lab = 0.5, pch = 12, 
 col = "pink", xlab = "variation price", ylab = "percent", 
 axes = TRUE, lablong = "x", lablat = "y")
}
}

\keyword{spatial}
\keyword{univar}
