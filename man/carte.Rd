\name{carte}
\alias{carte}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drawing a map}
\description{
The function \code{carte()} draws a map with sites of coordinates (\code{long},\code{lat}) and
represents sites which have been selected in \code{obs}.  This function is called by most of the functions of
GeoXp (this is not an interactive function).
}

\usage{
carte(long, lat, obs, sf.obj = NULL, num, criteria = NULL, buble = FALSE, 
 cbuble = NULL, nointer = FALSE, carte = NULL, nocart = FALSE, label = "",
 cex.lab = NULL, symbol = 16, lablong = "", lablat = "", method = "", 
 W = NULL, couleurs = "blue", classe = NULL, legmap = NULL, 
 legends = list(FALSE, FALSE), labmod = "", axis = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{long}{a vector x of size n}
  \item{lat}{a vector y of size n}
  \item{obs}{a boolean vector of size n with TRUE on selected sites}
  \item{sf.obj}{an sf object}
  \item{num}{number of the device}
  \item{criteria}{a vector of size n of boolean with TRUE on specific sites (these for non interactive selection)}
  \item{buble}{a boolean with TRUE for drawing bubbles, FALSE otherwise}
  \item{cbuble}{vector of size n with size of each site depending on variable with which bubbles are constructed }
  \item{nointer}{a boolean with TRUE for drawing sites selected by criteria, FALSE otherwise}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{nocart}{a boolean with TRUE for drawing spatial contours, FALSE otherwise}
  \item{label}{vector of character of size n with name of each site}
  \item{cex.lab}{character size of label}
  \item{symbol}{a vector of symbol which must be equal to the number of group else all sites are printed in pch[1]}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
  \item{method}{Specification for some graphics such as `Neighbourplot1', `Neighbourplot2', `Cluster', `Angleplot', `Variocloud'}
  \item{W}{neighbourhood matrix of size \eqn{n \times n}{n x n}, necessary when method=`Neighbourplot'}
  \item{couleurs}{Colors of the sites}
  \item{classe}{vector of n whith class of sites}
  \item{legmap}{a list of (numeric,numeric,numeric,numeric,numeric,numeric,character) with the sizes and the corresponding values of the 3 bubbles represented in the legend and the variable name of variable choosed }
  \item{legends}{a list of c(boolean,boolean,c(numeric,numeric),c(numeric,numeric)) with TRUE for drawing legends of bubbles (1st argument) or factors (2nd argument) and the coordinates of the location of the upper left corner of the legend box}
  \item{labmod}{Name of factors}
  \item{axis}{a boolean with True for drawing axes on the map}
}
\details{ Options permit to represent sites with different colors or symbols,  
with different sizes of circle (`bubble'), using a spatial contour (`carte'), printing names of sites (`label'),etc. 
Possibility to draw a legend giving the sizes of bubbles and colors of levels.
} 
\value{No values, only drawing of a map}
\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}
\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial} 
\keyword{utilities}
