\name{moranplotmap}
\alias{moranplotmap}
\title{Moran scatterplot and map}
\description{
The function \code{moranplotmap()} draws a moran plot, used to detect spatial
autocorrelation in the variable var. On the x-axis, is represented \eqn{x-\bar{x}}{x-bar(x)}
and on the y-axis \eqn{W(x-\bar{x})}{W(x-bar(x))}, where W is the spatial weight matrix.
It also calcultes Moran's I statistic (see \code{nonnormoran}) and
give a p-value associated to the autocorrelation test (gaussian version and permutation version). }

\usage{
moranplotmap(sf.obj, name.var, listw.obj, flower = FALSE, locmoran = FALSE, 
 names.arg = c("H.-H.", "L.-H.", "L.-L.", "H.-L."), 
 criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16, 
 col = "lightblue3", xlab = "", ylab = "", axes = FALSE, lablong = "", lablat = "")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{listw.obj}{object of class listw}
  \item{flower}{if TRUE, link neighbouring sites }
  \item{locmoran}{if TRUE, print local Moran's I statistic on the Moran plot}
  \item{names.arg}{names of the quadrant of the Moran plot}
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
For the permutation test, for each drawing the values of the variable \code{var} are
randomly assigned to the sites. We then calculate Moran'I statistic associated to each drawing and we give the frequency of drawings when Moran'I statistic is lower or equal to the
observed Moran'I statistic. Moreover, the function gives the opportunity to link neighbouring sites and gives Local Moran's I statistic.
For a site i:
\deqn{I_i=(x_i-\bar{x})\sum_{j}W_{ij}(x_j-\bar{x})}{I_i=(x_i-bar(x))sum_jW_ij(x_j-bar(x))}
with j not equal to i.}

\note{
In the case of the spatial weigth matrix is not normalized, the Moran'I statistic is not
equal to \eqn{\beta}{beta} used in regression line for model \eqn{W(X-\bar{X})=\beta(X-\bar{X})+u}{W(X-bar(X))=beta(X-bar(X))+u}. That is why
the regression line is only drawn in the case of W normalized.}

\value{
In the case where user click on \code{save results} button,
a list is created as a global variable in \code{last.select} object. \$obs, a vector of integer,
corresponds to the number of spatial units selected just before leaving the Tk window, \$MORAN,
a numeric, corresponds to the value of the Moran'I statistic.
}
\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer. \cr \cr

Jim Lesage, ``Spatial Econometrics Toolbox'', \url{http://www.spatial-econometrics.com/}}

\author{Aragon Y., Thomas-Agnan C., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}

\examples{
##
require(sf)
# data baltimore
if (require(spData, quietly = TRUE)) {
 data(baltimore, package = "spData")
}

# SpatialPoints
baltimore_sf <- st_as_sf(baltimore, coords = c("X", "Y"))
# Spatial Weight Matrix
W.nb <- spdep::knn2nb(spdep::knearneigh(baltimore_sf, k = 4))
# We choose a row standardized spatial weight matrix :
W.listw <- spdep::nb2listw(W.nb, style = "W")

# moranplotmap with some options
moranplotmap(baltimore_sf, "PRICE", W.listw, 
  flower = TRUE, locmoran = TRUE, criteria = (baltimore.spdf$AC == 1),
  identify = "STATION")

# comparison with the moran.test function
spdep::moran.test(baltimore_sf$PRICE, W.listw)
 
\dontrun{ 
# data columbus
require("sf")
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1])
col.gal.nb <- spdep::read.gal(system.file("weights/columbus.gal", package="spData")[1])

# use of moranplotmap with spatial weight matrix col.gal.nb :
# 1. row-standardized
moranplotmap(columbus, "HOVAL", nb2listw(col.gal.nb, style = "W"))

# 2. basic binary
moranplotmap(columbus,"HOVAL",nb2listw(col.gal.nb,style="B"))

# 3. globally standardized
moranplotmap(columbus, "HOVAL", nb2listw(col.gal.nb, style="C"))
}
}


\seealso{\code{\link{neighbourmap}},\code{\link{makeneighborsw}},\code{\link{normw}},\code{\link{nonormmoran}} }
