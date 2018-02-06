\name{polylist2list}
\alias{polylist2list}

\title{Extract from a polylist object the vertices of the polygons}

\description{
The function \code{polylist2list()} extracts the vertices of the polygons from a \code{polylist} object, seperating polygons from each other by 3 rows of NaN.
}
\usage{
polylist2list(data)
}

\arguments{
  \item{data}{A \code{polylist} object}
}

\details{
The user can then represent the coordinates of sites of a \code{polylist} object
as background map, using the option \code{carte} included in all interactive functions of GeoXp.
}

\value{
It returns a matrix of numeric values with 2 columns (x and y coordinates
of the vertices of the polygons) where polygons are seperated from each other by 3 rows of NaN.
}

\references{
Bivand R. et al. (2009), \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{T. Laurent}

\seealso{\code{\link{spdf2list}}}



\examples{
data(eire)
eire.contours<-polylist2list(eire.polys.utm)
}

\keyword{manip}
