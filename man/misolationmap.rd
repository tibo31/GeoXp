\name{misolationmap}
\alias{misolationmap}
\title{Interactive multivariate isolation plot and map}
\description{
The function \code{misolationmap} draws a scatterplot with the pairwise Mahalanobis distances calculated using
variables \code{names.var} between the observations and their neighbors on the y-axis and the "degree of isolation" of the observations on
the x-axis and a map
}
\usage{
misolationmap(sp.obj, nb.obj, names.var, propneighb = 0.4, chisqqu = 0.975,
 names.attr = names(sp.obj), criteria = NULL, carte = NULL, identify = FALSE, 
 cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "degree of isolation", 
 ylab = "Pairwise Mahalanobis distances", axes = FALSE, 
 lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{nb.obj}{object of class nb}
  \item{names.var}{a vector of character; attribute names or column numbers in attribute table}
  \item{propneighb}{proportion of neighbors included in ellipsoid}
  \item{chisqqu}{value of alpha for the definition of global outliers}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of size n of boolean which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{color of the points on the cloud map}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}

\details{
The pairwise Mahalanobis distances are calculated using the robust Minimum Covariance Determinant (MCD)
estimator associated with \eqn{75\%}{75\%} of observations (function \code{covMcd} in the robustbase package)
calculated on the variables \code{names.var}.
For each observation, the degree of isolation is a chi-square quantile of the conditional distribution of the
pairwise Mahalanobis distances associated with the ellipsoid containing the proportion \code{propneighb}
of neighbors. The parameter \code{propneighb} gives the proportion of neighbors that is expected to be quite
similar to the observation in order to conclude that the observation is not a local outlier.
Under independence and normality conditions, the user can expect a degree of isolation close by the parameter
\code{propneighb} (vertical line on the scatterplot). An observation with a high degree of isolation is suspected
to be a local outlier. Users have also the possibility to plot bubbles on the map which size depends on the
robust Mahalanobis distance of each observation to the center of the distribution (function \code{arw} in the
package \code{mvoutlier}).
}

\value{
In the case where user click on \code{save results} button,
a matrix of integer is created as a global variable in \code{last.select} object.
It corresponds to the numbers of spatial unit corresponding to couple of sites selected
just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Fizmoser P., Thomas-Agnan C., Ruiz-Gazen A., Laurent T.,}

\examples{
## data radarImage
require("robustbase")
data(radarImage)

# creation of a SpatialClass object
radarImage.sp <- SpatialPoints(radarImage[1300:1573, c("X.coord", "Y.coord")])
radarImage.spdf <- SpatialPointsDataFrame(radarImage.sp, radarImage[1300:1573,])

# creation of a spatial weight matrix nb
radarImage.nb <- dnearneigh(radarImage.sp, 0, 1.5)

# example of use of misolationmap
# The statistics are calculated by taking into account variables
# Ag,As,Bi,Cd,Co,Cu,Ni
misolationmap(radarImage.spdf, radarImage.nb, names.var = c("Band.1","Band.2","Band.3"),
 propneighb = 0.30, chisqqu = 0.95, identify = TRUE, cex.lab = 0.5)

}

\keyword{spatial}
\keyword{multivariate}
\seealso{\code{\link{mvariocloudmap}}}
