\name{mvariocloudmap}
\alias{mvariocloudmap}
\title{Interactive multivariate variocloud and map}
\description{
The function \code{mvariocloudmap()} draws a scatterplot of pairwise Mahalanobis
distances and spatial distances with a map. It is a multivariate version of the
variocloud. The number of couples of sites plotted can be reduced by considering
couples above a quantile regression curve.
}
\usage{
mvariocloudmap(sp.obj, nb.obj, names.var, quantiles=TRUE, 
names.attr=names(sp.obj), criteria = NULL, carte = NULL, identify = FALSE, cex.lab = 0.8,
 pch = 16, col = "lightblue3", xlab = "Pairwise spatial distances", 
 ylab = "Pairwise Mahalanobis distances", axes = FALSE, lablong = "", lablat = "")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{nb.obj}{object of class nb}
  \item{names.var}{a vector of character; attribute names or column numbers in attribute table}
  \item{quantiles}{a boolean to represent the Additive Quantile Regression Smoothing}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean of size the number of spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
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
The pairwise Mahalanobis distances are calculated using the Minimum Covariance Determinant (MCD)
estimator associated with \eqn{75\%}{75\%} of observations (function \code{covMcd} in the \code{robustbase}
package). Users have the possibility to select some couples of sites on the scatterplot that are also highlightened
on the map. Selection of observations on the map is also possible and leads to the selection of all the
couples which contain the selected observations on the scatterplot.
}

\value{
In the case where user click on \code{save results} button,
a matrix of integer is created as a global variable in \code{last.select} object.
It corresponds to the numbers of spatial unit corresponding to couple of sites selected
just before leaving the Tk window.
}


\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Fizmoser P., Thomas-Agnan C., Ruiz-Gazen A., Laurent T.}

\examples{
## data meuse
data(meuse)

# transformation of explanatory variables
meuse[, 3:7] <- log(1 + meuse[, 3:7])

# creation of a Spatial Points object
meuse.sp <- SpatialPoints(cbind(meuse$x, meuse$y))

# creation of a SpatialPointsDataFrame
meuse.spdf <- SpatialPointsDataFrame(meuse.sp, meuse)

# for the spatial contours
data(meuse.riv)

# Spatial Weight matrix based on the 7th nearest neighbours
meuse.knn <- knearneigh(meuse.sp, k = 7)
meuse.nb <- knn2nb(meuse.knn)

# example of use of mvariocloudmap. The statistic are calculated by taking
# into account variables cadmium,copper,lead,zinc,elev
\donttest{
mvariocloudmap(meuse.spdf, meuse.nb, c("cadmium", "copper", "lead", "zinc", "elev"), 
 quantiles = 0.95, carte = meuse.riv[-c(1:20, 73:98, 156:176), ], identify = TRUE, 
 criteria = (meuse.spdf$lime == 1))
 }
}

\keyword{spatial}
\keyword{multivariate}
\seealso{\code{\link{misolationmap}}}
