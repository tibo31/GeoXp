\name{clustermap}
\alias{clustermap}
\title{Classification of dataset using kmeans or hclust algorithm and representation of clusters on a map.}
\description{
The function \code{clustermap()} performs a classification of the sites from the variables called
in \code{names.var} and computes a bar plot of the clusters calculated.
Classification methods come from \code{hclust()} (hierarchical cluster analysis) and \code{kmeans()}
(k-means clustering) and number of class is chosen with \code{clustnum}.
}

\usage{
clustermap(sf.obj, names.var, clustnum, method = c("kmeans", "hclust"), 
  type = NULL, centers = NULL, scale = FALSE, names.arg = "", 
  criteria = NULL, carte = NULL, identify = NULL, 
  cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "Cluster", 
  ylab = "Number", axes = FALSE, lablong = "", lablat = "")
}

\arguments{
  \item{sf.obj}{object of class sf}
  \item{names.var}{a vector of character; attribute names or column numbers in attribute table}
  \item{clustnum}{integer, number of clusters}
  \item{method}{two methods : "kmeans"" by default or "hclust"}
  \item{type}{If method="hclust", type="complete" by default  (the possibilities are given in help(hclust)
  as "ward", "single", etc). If method="kmeans", type="Hartigan-Wong" by default (the possibilities are given in help(kmeans)
  as "Forgy", etc)}
  \item{centers}{If method='kmeans', user can give a matrix with initial cluster centers.}
  \item{scale}{If scale=TRUE, the dataset is reducted.}
  \item{names.arg}{a vector of character, names of cluster}
  \item{criteria}{a vector of boolean of size the number of spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
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
The two windows are interactive : the sites selected by a bar chosen on the bar plot 
are represented on the map in red and the values of sites selected on the map by `points' or `polygon' are 
represented in red on the bar plot. The dendogram is also drawn for 'hclust' method.  
In option, possibility to choose the classification method. 
}
\value{
In the case where user click on \code{save results} button,
a list is created as a global variable in \code{last.select} object. \code{obs}, a vector of integer,
corresponds to the number of spatial units selected just before leaving the Tk window, \code{vectclass},
vector of integer, corresponds to the number of cluster attributed to each spatial unit.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Murtagh, F (1985). \emph{Multidimensional Clustering Algorithms}. \cr \cr

Hartigan, J. A. and Wong, M. A. (1979). A K-means clustering algorithm. Applied Statistics 28, 100-108 \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\note{To use the functions \code{hclust} and \code{kmeans}, we take many arguments by default.
If the user would like to modify these arguments, he should call these functions first and
then use the function \code{barmap} to visualize the calculated clusters. }

\seealso{\code{\link{barmap}}, \code{\link{pcamap}}}

\examples{
#####
# data columbus
require("sf")
columbus <- st_read(system.file("shapes/columbus.shp", package="spData")[1])

# a basic example using the kmeans method
clustermap(columbus, c("HOVAL", "INC", "CRIME", "OPEN", "PLUMB", "DISCBD"), 3, 
  criteria = (columbus@data$CP == 1), identify = TRUE, cex.lab = 0.7)

\dontrun{ 
# example using the hclust method
clustermap(columbus, c(7:12), 3, method = "hclust", 
  criteria = (columbus$CP == 1), col = colors()[20:22], identify = "POLYID", 
  cex.lab = 0.7, names.arg = c("Group 1", "Group 2", "Group 3"), xlab = "Cluster")
}
}

\keyword{spatial}
\keyword{multivariate}
