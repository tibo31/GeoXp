\name{selectmap}
\alias{selectmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of a point or polygon on a scatterplot}

\description{
The function \code{selectmap()} is used to select a point or a polygon on the map or
on a scatterplot. Called by any function which draws a scatterplot.
}

\usage{
selectmap(var1, var2, obs, Xpoly, Ypoly, method = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var1}{a vector x of size n}
  \item{var2}{a vector y of size n }
  \item{obs}{a boolean vector of size n : TRUE if a site is already selectioned, FALSE otherwise.}
  \item{Xpoly}{X-coordinates of the vertices of selected polygon}
  \item{Ypoly}{Y-coordinates of the vertices of selected polygon}
  \item{method}{`point' if the selected area is a point; `poly' if the selected area is a polygon}
}

\details{This function is called by all the functions which draw a scatterplot such as \code{scattermap},
\code{moranplotmap},... This is not an interactive function.}

\value{
A vector of boolean of size n. TRUE if a site has been selected, FALSE otherwise.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y.,  Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}
\keyword{utilities} 