\name{selectstat}
\alias{selectstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of values on a graphic}

\description{The function \code{selectstat()} is used to select sites on a graphic such
as Histogram, Barplot, Boxplot, Polyboxplot, Scatterplot, Densityplot, Neighbourplot, Lorentz,
Moran, Angleplot, Variocloud, Cluster, or PCA. Most of the GeoXp functions call selectstat.r.
}
\usage{
selectstat(var1, var2, obs, Xpoly, Ypoly, method,
nbcol, W, F, long, lat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var1}{1st variable of size n }
  \item{var2}{2nd variable of size n }
  \item{obs}{a boolean vector of size n : TRUE if a site is already selectioned, FALSE otherwise. }
  \item{Xpoly}{X-coordinates of the vertices of selected polygon}
  \item{Ypoly}{Y-coordinates of the vertices of selected polygon}
  \item{method}{name of the graph wich must be drawn among Histogram, Barplot, Boxplot, Polyboxplot, Densityplot, Neighbourplot, Lorentz, Anglepoint, Variopoint}
  \item{nbcol}{nbcol : number of cells for the histogram (10 by default)}
  \item{W}{A matrix of size \eqn{n \times n}{n x n} (for spatial weight)}
  \item{F}{a vector of numeric}
  \item{long}{a vector of x-axis of size n}
  \item{lat}{a vector of y-axis of size n}
}
\details{This function is not an interactive function.}
\value{Return a vector of boolean of size n with TRUE if sites have been selected on a graph}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}
\keyword{utilities}