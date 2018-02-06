\name{histobarmap}
\alias{histobarmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Histogram, barplot and map}

\description{
The function \code{histobarmap()} draws a bar plot (vertical bar) of the given variable
\code{names.var[1]}, a histogram of the given variable \code{names.var[2]} and a map with
sites of coordinates \code{coordinates(sp.obj)}.}

\usage{
histobarmap(sp.obj, names.var, nbcol = 10, type = "count",
 names.arg = "", names.attr = names(sp.obj), criteria = NULL, carte = NULL, 
 identify = FALSE, cex.lab = 0.8, pch = 16, col = "lightblue3", 
 xlab = c("barplot", "histogram"), ylab = rep("count", 2),
 axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{names.var}{a vector of 2 characters; first character corresponds to the name of a factor and second character corresponds to a numeric value}
  \item{nbcol}{number of cells for histogram (10 by default)}
  \item{type}{Character string indicating type of histogram/barplot to be drawn. "percent" and "count" give relative frequency and frequency histogram/barplot.}
  \item{names.arg}{a vector of level names (for factor)}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean of size the number os spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
  \item{cex.lab}{character size of label}
  \item{pch}{a vector of symbol which must be equal to the number of group else all sites are printed in pch[1]}
  \item{col}{a vector of colors which must be equal to the number of group else all sites and all box are printed in col[1]}
  \item{xlab}{a vector of names for x-title of the barplot and histogram}
  \item{ylab}{a vector of names for y-title of the barplot and histogram}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}
\details{
Each site is associated to a factor of \code{names.var[1]} and to a value of
\code{names.var[2]}. There is interactivity between the three windows created:
the sites selected by a bar on the bar plot or on the histogram are represented on the
map in red and the value and factor of sites selected on the map are represented in red
on the bar plot and on the histogram.
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan C., Aragon Y.,  Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{dblehistomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}},  \code{\link{dbledensitymap}}}



\examples{
###
# Data Colombus
require("rgdal")
columbus <- readOGR(system.file("shapes/columbus.shp", package="spData")[1])

# an example of use
histobarmap(columbus, c("CP", "HOVAL"), nbcol = 8, type = "percent", 
 names.arg = c("A", "B"), xlab = c("CP", "Hoval"), ylab = rep("percent", 2))

\dontrun{ 
# data eire
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])

# example of use
histobarmap(eire, c("pale", "A"),names.arg = c("Ouside Pale", "Pale"),
 xlab = c("Appartenance to the region of Pale", "Average number of people with blood A"),
 col = colors()[101:102], identify = TRUE)
}
}

\keyword{spatial}
\keyword{multivariate} 