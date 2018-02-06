\name{makeneighborsw}
\alias{makeneighborsw}
\title{Spatial weight matrix}
\description{
The function \code{makeneighborsw()} create a spatial weight matrix based on a given number of nearest
neighbors (option "neighbor" by default), based on a threshold distance (option method="distance")
or both these 2 methods.}

\usage{
makeneighborsw(coords,method="neighbor",m=1,d,cum=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coords}{a matrix of spatial coordinates}
  \item{method}{"neighbor" by default, "distance" or "both"}  
  \item{m}{number of nearest neighbors}
  \item{d}{threshold point}
  \item{cum}{if cum=TRUE, W is the sum of spatial weight matrix based on k nearest neighbours (for \eqn{k \leq m}{k < m};
   if FALSE W is the spatial weight matrix based only on \eqn{m^{th}}{m^th} nearest neighbours}
}
\details{
In the case of method="neighbor", for each site, we order the other sites by their distance from this
site. If cum=TRUE, for i, if j is among the \eqn{m^{th}}{m^th}
nearest sites, then :     
\deqn{W_{ij}=1}{W_ij=1}
else
\deqn{W_{ij}=0}{W_ij=0}
If cum=FALSE, for \deqn{s_i}{s_i}, if \deqn{s_j}{s_j} is the \eqn{m^{th}}{m^th}
nearest site, then :     
\deqn{W_{ij}=1}{W_ij=1}
else
\deqn{W_{ij}=0}{W_ij=0}
In case of ties, the nearest neighbour is randomly chosen. \cr \cr
In the case of method="distance", if site i is seperated from j by a distance lower
or equal to a given threshold : 
\deqn{W_{ij}=1}{W_ij=1}
else
\deqn{W_{ij}=0}{W_ij=0}
\cr 
In the case of method="both" W must verify the two first conditions.}

\note{
This function is not optimised for large dataset. User could find similar functions in
the package \code{spdep} (\code{dnearneigh} and \code{knearneigh}). However, these functions
don't offer the possibility to use the two criteria in the same time.  Moreover, an
inconvenient of \code{makeneighborsw} is that the result is included in a matrix object
whereas most of functions of GeoXp use the \code{nb} structure for spatial weight matrix.
An issue is to use the \code{mat2listw} function and then selecting the \code{nb} part, like the
in the examples. }

\value{
A spatial weight matrix of size \eqn{n \times n}{n x n}}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Aragon Y., Thomas-Agnan C., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}

\seealso{\code{\link{moranplotmap}},\code{\link{normw}}}

\examples{
##
# data auckland
data(auckland)
coords <- cbind(auckland$Easting[1:10],auckland$Northing[1:10])

# matrix based on 5 nearest neighbors
W<-makeneighborsw(coords,method="neighbor",m=3)

# matrix based on a threshold distance
W1<-makeneighborsw(coords,method="distance",d=20)

# matrix based on the two methods
W2<-makeneighborsw(coords,method="both",m=3,d=20)

# representation of the 3 spatial weight matrices
op<-par(mfrow=c(2,2))
plot(mat2listw(W),coords,col="lightblue1",main="neighbor")
plot(mat2listw(W1),coords,col="lightblue2",main="distance")
plot(mat2listw(W2),coords,col="lightblue3",main="both")
par(op)
} 