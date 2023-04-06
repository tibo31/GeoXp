\name{variocloudmap}
\alias{variocloudmap}
\title{Interactive variocloud and map}

\description{
The function \code{variocloudmap()} draws a semi-variocloud (directional or omnidirectional) and a map.
It is used to detect spatial autocorrelation. Possibility to draw the empirical semi-variogram
and a robust empirical semi-variogram.}

\usage{
variocloudmap(sf.obj, name.var, bin = NULL, quantiles = TRUE,
  criteria = NULL, carte = NULL, identify = NULL, 
  cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", axes = FALSE,
  lablong = "", lablat = "", xlim = NULL, ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{bin}{a vector of numeric values where empirical variogram will be evaluated}
  \item{quantiles}{a boolean to represent the Additive Quantile Regression Smoothing}
  \item{criteria}{a vector of boolean of size the number of Spatial Units, which permit to represent preselected sites with a cross, using the tcltk window}
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
  \item{xlim}{the x limits of the plot}
  \item{ylim}{the y limits of the plot}
}

\details{
For some couple of sites \eqn{(s_i,s_j)}{(s_i,s_j)}, the graph represents on the y-axis the semi squared difference
between \eqn{var_i}{var_i} and \eqn{var_j}{var_j} :
\deqn{\gamma_{ij}=\frac{1}{2}(var_i-var_j)^2}{gamma_ij=0.5(var_i-var_j)^2}
and on the x-absis the distance \eqn{h_{ij}}{h_(ij)} between \eqn{s_i}{s_i} and \eqn{s_j}{s_j}.
The semi Empirical variogram has been calculated as :
\deqn{\gamma(h)=\frac{1}{2|N(h)|}\sum_{N(h)}(Z(s_i)-Z(s_j))^2}{gamma(h)=0.5/|N(h)|sum_(N(h))(Z(s_i)-Z(s_j))^2}
where
\deqn{N(h)=\{(s_i,s_j):s_i-s_j=h;i,j=1,...,n\}}{N(h)=\{(s_i,s_j):s_i-s_j=h;i,j=1,...,n\}}
and the robust version :
\deqn{\gamma(h)=\frac{1}{2(0.457+\frac{0.494}{|N(h)|})}(\frac{1}{|N(h)|}\sum_{N(h)}|Z(s_i)-Z(s_j)|^{1/2})^4}{gamma(h)=frac(1)(2(0.457+frac(0.494)(|N(h)|)))(frac(1)(|N(h)|)sum_(N(h))|Z(s_i)-Z(s_j)|^(1/2))^4}
The number N of points to evaluate the empirical variogram and the distance \eqn{\epsilon}{epsilon} between
points are set as follows :
\deqn{N=\frac{1}{max(30/n^2,0.08,d/D)}}{N=frac(1)(max(30/n^2,0.08,d/D))}
and :
\deqn{\epsilon=\frac{D}{N}}{epsilon=frac(D)(N)}
with : \deqn{D=max(h_{ij})-min(h_{ij})}{D=max(h_ij)-min(h_ij)}
and : \deqn{d=max(h_{ij}^{(l)}-h_{ij}^{(l+1)}),}{d=max(h_ij^(l)-h_ij^(l+1)),}
where \eqn{h^{(l)}}{h^(l)} is the vector of sorted distances.
In options, possibility to represent a regression quantile smoothing spline \eqn{g_{\alpha}}{g_alpha} (in that
case the points below this quantile curve are not drawn). }

\value{
In the case where user click on \code{save results} button,
a matrix of integer is created as a global variable in \code{last.select} object.
It corresponds to the numbers of spatial unit corresponding to couple of sites selected
just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Cressie N. and Hawkins D. (1980), \emph{Robust estimation of the variogram}, in Journal of the international association for mathematical geology, 13, 115-125.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}
\examples{
#####
require(sf)
# Data Meuse
if (require(sp, quietly = TRUE)) {
 data(meuse, package = "sp")
 meuse_sf <- st_as_sf(meuse, coords = c("x", "y"), crs = 28992, agr = "constant")
 data(meuse.riv, package = "sp")
 meuse.sr <- st_as_sf(SpatialPolygons(list(Polygons(list(Polygon(meuse.riv)), 
 "meuse.riv"))), crs = 28992)
}

# example of use of variocloudmap
\donttest{
variocloudmap(meuse_sf, "zinc", quantiles = TRUE, bin = seq(0, 2000, 100),
  xlim = c(0, 2000), ylim = c(0, 500000), pch = 2, carte = meuse.sr,
  criteria = (meuse_sf$lime == 1))
 }
}

\seealso{\code{\link{angleplotmap}}, \code{\link{driftmap}}}