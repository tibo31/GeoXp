\name{angleplotmap}
\alias{angleplotmap}
\title{Detection of an eventual directional trend}
\description{
The function \code{angleplotmap()} is used to detect an eventual directional trend associated
to variable \code{name.var}. It represents the absolute difference between the value of \code{name.var}
at two sites as a function of the angle between vector \eqn{\overrightarrow{s_is_j}}{vector(s_is_j)}
and the x-axis.
}
\usage{
angleplotmap(sf.obj, name.var, quantiles = TRUE, 
  criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16, 
  col = "lightblue3", xlab = "angle", ylab = "absolute magnitude", axes = FALSE,
  lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{quantiles}{a boolean to represent the Additive Quantile Regression Smoothing}
  \item{criteria}{a vector of boolean which permit to represent preselected sites with a cross, using the tcltk window (must be equal to the number of spatial units)}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{"lightblue3" by default, color of bars on the cloud map}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
For each couple of sites \eqn{(s_i,s_j)}{(s_i,s_j)}, the graphic represents on the y-axis the absolute difference
between \eqn{var_i}{var_i} and \eqn{var_j}{var_j} :
\deqn{D_{ij}=|var_i-var_j|}{D_ij=abs(var_i-var_j)} \cr and on the x-axis the angle \eqn{\theta_{ij}}{theta_ij} between
\eqn{\overrightarrow{s_is_j}}{vector(s_is_j)} and the x-axis.
Possibility to represent a smoothing spline regression quantile \eqn{g_\alpha}{g_alpha}. For \eqn{0<\alpha<1}{0<alpha<1},
\deqn{Pr[D_{ij}<g_\alpha(\theta{ij})]=\alpha}{Pr[D_ij<g_alpha(theta_ij)]=alpha}
If that case, only the pair of sites \eqn{(s_i,s_j)}{(s_i,s_j)} verifying :
\deqn{D_{ij}>g_{max(\alpha)}(\theta{ij})}{D_ij>g_max(alpha)(theta_ij)} are represented.
}

\value{
In the case where user click on \code{save results} button,
a matrix of integer is created as a global variable in \code{last.select} object.
It corresponds to the numbers of spatial unit corresponding to couple of sites selected
just before leaving the Tk window.
}


\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan Christine, Aragon Yves, Ruiz-Gazen Anne, Laurent Thibault, Robidou Lauriane}

\keyword{spatial}

\examples{
######
# data Conflicts Africa
if (require(spData, quietly = TRUE)) {
 data(afcon, package = "spData")
}
# afcon is a data.frame object. We have to create
# a Spatial object, by using first the longitude and latitude
# to create Spatial Points object ...
require(sf)
afcon.sf <- st_as_sf(afcon, coords = c("x", "y"), crs = 4326)

# optional : we add some contours that don't correspond to the spatial unit
# but are nice for mapping
africa <- st_read(system.file("shapes/Africa.shp", package = "GeoXp")[1])

# A basic call of histomap function
angleplotmap(afcon.sf, "totcon", carte = africa,
  identify = "name", cex.lab = 0.6)
 
\dontrun{ 
if (require(sp, quietly = TRUE)) {
 data(meuse, package = "sp")
 meuse_sf <- sf::st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
 data(meuse.riv, package = "sp")
 meuse.sr <- sf::st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))),
     crs = 28992)
}

angleplotmap(meuse_sf, "copper", col = "green", quantiles = TRUE, 
  cex.lab = 0.7, xlab = "Concentration in plomb (in ppm)", pch = 7, 
  carte = meuse.sr)
}
}

\seealso{\code{\link{variocloudmap}}, \code{\link{driftmap}}}