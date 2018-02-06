\name{normw}
\alias{normw}

\title{Row-normalize a spatial weight matrix}
\description{
The function \code{normw()} row-normalizes a spatial weight matrix}

\usage{
normw(w)
}

\arguments{
  \item{w}{A matrix of size \eqn{n \times n}{n x n}}
}
\details{
 \deqn{W_{ij}=\frac{W_{ij}}{\sum_kW_{ik}}}{W_ij=frac(W_(ij))(sum_kW_(ik))}
}

\value{
A matrix of size \eqn{n \times n}{n x n}
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Aragon Y., Thomas-Agnan C., Ruiz-Gazen A., Robidou L.}

\keyword{spatial}

\seealso{\code{\link{neighbourmap}}, \code{\link{moranplotmap}},\code{\link{makeneighborsw}}}

\examples{
###
# data auckland
data(auckland)
x.ext <- auckland$Easting[1:10]
y.ext <- auckland$Northing[1:10] 

# matrix based on 3 nearest neighbors
W1<-makeneighborsw(cbind(x.ext,y.ext),method="both",m=3,d=20)
W2<-normw(W1)

apply(W1[1:10,],1,sum)
apply(W2[1:10,],1,sum)
} 