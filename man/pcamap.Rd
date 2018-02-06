\name{pcamap}
\alias{pcamap}
\title{Generalized Principal Component Analysis and map}
\description{
The function \code{pcamap()} draws the plots summarizing a generalized Principal Component Analysis (PCA),
made with \code{genpca}. It draws the scatterplot of the individuals projected on a chosen principal
component plane (with their percentage of inertia), together with the scatterplot of the variables
projected into the same plane with the quality of representation in order to interpret the principal
component axes. The individuals scatterplot interacts with the map.
}
\usage{
pcamap(sp.obj, names.var, direct = c(1, 2), weight = rep(1/nrow(sp.obj),
 length = nrow(sp.obj)), metric = diag(length(names.var)), 
 center = NULL, reduce = TRUE, qualproj = FALSE, names.attr = names(sp.obj), 
 criteria = NULL, carte = NULL, identify = FALSE, cex.lab = 0.8, pch = 16, 
 col = "lightblue3", xlab = paste(direct[1]), ylab = paste(direct[2]), 
 axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{names.var}{a vector of character; attribute names or column numbers in attribute table}
  \item{direct}{Two-dimension vector containing the numbers of principal axes to plot}
  \item{weight}{vector of size n of weight (by default : weight=t(1/n,...,1/n))}
  \item{metric}{ matrix \eqn{p \times p}{p x p} (by default : metric=Identity matrix)}
  \item{center}{A vector of size p of the gravity center which by default is equal to \eqn{\bar X}{bar(X)}}
  \item{reduce}{if TRUE, reduced PCA}
  \item{qualproj}{if TRUE, print the quality of representation of individuals }
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean  of size the number of Spatial unit, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
  \item{cex.lab}{character size of label}
  \item{pch}{a vector of symbol which must be equal to the number of group else all sites are printed in pch[1]}
  \item{col}{a vector of colors which must be equal to the number of group else all sites and all bars are printed in col[1]}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
Let \deqn{D=diag(\lambda_1,...,\lambda_p)}{D=diag(lambda_1,...,lambda_p)}
\deqn{1_p=(1,...,1)'}{1_p=(1,...,1)'}
Let the coordinates of individuals in the principals components \deqn{CC=(C_1',...,C_n')'}{CC=(C_1',...,C_n')'}
with \eqn{C_i=(C_i^1,...,C_i^p)}{C_i=(C_i^1,...,C_i^p)} \cr
Let the coordinates of variables in the principals components \deqn{CC=(V_1',...,V_p')'}{CC=(V_1',...,V_p')'}
with \eqn{V_i=(V_i^1,...,V_i^p)}{V_i=(V_i^1,...,V_i^p)} \cr
Part of inertia : \deqn{(\frac{\lambda_1}{\sum_i\lambda_i},...,\frac{\lambda_p}{\sum_i\lambda_i})'}{(frac(lambda_1)(sum_ilambda_i),...,frac(lambda_p)(sum_ilambda_i))'}
Quality of representation of individual k projected on plane (i,j):
\deqn{Qu=\sqrt{\frac{(C_k^i)^2+(C_k^j)^2}{\sum_l(C_k^l)^2}}}{Qu=sqrt(frac((C_k^i)^2+(C_k^j)^2)(sum_l(C_k^l)^2)}
Quality of representation of variable k projected on plane (i,j):
\deqn{VQu=\sqrt{\frac{(V_k^i)^2+(V_k^j)^2}{\sum_l(V_k^l)^2}}}{VQu=sqrt(frac((V_k^i)^2+(V_k^j)^2)(sum_l(V_k^l)^2))}
}

\value{
In the case where user click on \code{save results} button,
a list is created as a global variable in \code{last.select} object. \code{obs},
corresponds to the number of spatial units selected just before leaving the Tk window,
\code{inertia} vector of size p with percent of inertia of each component, \code{casecoord}
matrix \eqn{n \times p}{n x p} of individuals, \code{varcoord} matrix \eqn{n \times p}{n x p}
of principal components.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Caussinus H., Fekri M., Hakam S., Ruiz-Gazen A. (2003) , \emph{A monitoring display of Multivariate Outliers}, Computational Statistics and Data Analysis, vol. 44, 1-2, 237-252.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\examples{
###
# Data Colombus
require("robustbase")
require("rgdal")
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])


# a basic PCA on 7 variables
pcamap(columbus, c(6:12), qualproj = TRUE, identify = TRUE)

###
# data boston
data(boston)

# SpatialPoints object
boston.sp <- SpatialPoints(cbind(boston.utm[, 1], boston.utm[, 2]))

# SpatialPointsDataFrame object
boston.spdf <- SpatialPointsDataFrame(boston.sp, boston.c)

# a basic PCA on 7 variables
pcamap(boston.spdf, c(7:8, 10, 12:15), identify = TRUE, cex.lab = 0.5)

# generalized PCA : user have to construct a new metric and a vector
# of gravity center, by using for exampe covMcd
cov.boston <- covMcd(boston.c[, c(7:8, 10, 12:15)], alpha = .75)
b.center <- cov.boston$center
b.cov <- cov.boston$cov

# example of use of pcamap
pcamap(boston.spdf, c(7:8, 10, 12:15), metric = b.cov, center = b.center,
 identify = TRUE, cex.lab = 0.5)
}

\keyword{spatial}
\keyword{multivariate}
\seealso{\code{\link{clustermap}},\code{\link{genpca}}}
