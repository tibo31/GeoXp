\name{driftmap}
\alias{driftmap}

\title{Interactive driftplot and map}
\description{
The function \code{driftmap()} creates two devices : a device with the map of coord
\code{coordinates(sf.obj)} which permits to make selection of sites and a device divided into
2 rows and 2 columns which contains : (cell 1) the selected sites
divided into m rows and q columns (m and q are selected with the tk window),
(cell 2) a scatter plot with \code{coordinates(sf.obj)[,2]}
in x-axis and the mean and median of \code{name.var} calculated for the m rows in y-axis,
(cell 3) a scatter plot with the mean and
median of \code{name.var} calculated for the q columns in x-axis and \code{coordinates(sp.obj)[,1]} in y-axis and (cell 4)
a legend indicating the direction of the North, the variable used, etc.}

\usage{
driftmap(sf.obj, name.var, interpol = TRUE, nuage = TRUE, lty = 1:2, 
  cex = 0.7, carte = NULL, identify = "", 
  cex.lab = 0.8, pch = rep(16,3), col = c("lightblue3", "black", "red"), 
  xlab = "", axes = FALSE)
}

\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{interpol}{if TRUE, the mean and median calculated are linearly interpoled}
  \item{nuage}{if TRUE, the values taken by var are also represented on right plot and on the left plot}
  \item{lty}{the line type for mean and median in the case of interpol=TRUE}
  \item{cex}{the amount by which plotting symbols on the cell 1 and in the case of nuage=TRUE cell 2 and cell 3}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{c(16,16,16) by default, list of symbols for (1) the values of var on the left and on the right plot (2) the mean points, (3) the median points}
  \item{col}{c("lightblue3", "black", "red") by default, list of colors for (1) the values of var on the left and on the right plot (2) the mean points, (3) the median points}
  \item{xlab}{name of var printed with the legend plot}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
}

\details{
Possibility to change the number of cells in the grid with the tk window, to interpolate the means and medians calculated (by default),
to work on a rotated map. At the begining, all sites have been selected to the map : users have to deselect sites that they
have choosen.}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.
}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{angleplotmap}},\code{\link{variocloudmap}},\code{\link{rotation}}}

\examples{
######
# data on price indices of real estate in France
data(immob, package = "GeoXp")

# immob is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
require(sf)
immob.sf <- st_as_sf(immob, coords = c("longitude", "latitude"))
# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
midiP <- st_read(system.file("shapes/region.shp", package="GeoXp")[1])

# a simple use of driftmap
driftmap(immob.sf, "prix.vente", carte = midiP)

# ... with options
driftmap(immob.sf, "prix.vente", carte = midiP, 
  interpol = FALSE, nuage = TRUE, lty = 3:4, identify = "Nom",
  cex.lab = 0.3, xlab = "Price of sell", axes = TRUE)

\dontrun{ 
# data eire
require(sf)
eire <- st_read(system.file("shapes/eire.shp", package="spData")[1])

# a basic usage ...
driftmap(eire, "pale", identify = "names", nuage = TRUE)
}
}

\keyword{spatial}
