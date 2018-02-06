\name{barmap}
\alias{barmap}
\title{Interactive Bar plot and map}
\description{
The function \code{barmap()} draws a bar plot (vertical bar) of the given factor variable
\code{name.var} and a map with sites of coordinates \code{coordinates(sp.obj)}.}

\usage{
barmap(sp.obj, name.var, type = c("count","percent"), names.arg = "", 
 names.attr = names(sp.obj), criteria = NULL, carte = NULL, identify = FALSE,
 cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", 
 axes = FALSE, lablong = "", lablat = "")
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sp.obj}{object of class extending Spatial-class}
  \item{name.var}{a character; attribute name or column number in attribute table}
  \item{type}{Character string indicating type of histogram to be drawn. "percent" and "count" give relative frequency and frequency histograms.}
  \item{names.arg}{a vector of level names of name.var}
  \item{names.attr}{names to use in panel (if different from the names of variable used in sp.obj)}
  \item{criteria}{a vector of boolean of size the number os spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not FALSE, identify plotted objects (currently only working for points plots). Labels for identification are the row.names of the attribute table row.names(as.data.frame(sp.obj)).}
  \item{cex.lab}{character size of label}
  \item{pch}{a vector of symbol which must be equal to the number of level else all sites are printed whith pch[1]}
  \item{col}{a vector of colors which must be equal to the number of level else all sites and all bars are printed with col[1]}
  \item{xlab}{a title for the graphic x-axis}
  \item{ylab}{a title for the graphic y-axis}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}  
}
\details{The selection of a bar on the bar plot results in the corresponding 
sites coloured on the map with the corresponding colour observed on the bar.
Reversely, the selection of sites on the map by "points" or "polygon"
results in the drawing of the sub-barplot in red.
}

\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{histomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}

\examples{
######
# data eire
require("rgdal")
eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])

# a basic usage ...
barmap(eire, "pale")

# ... with all options
barmap(eire, 3, type = "percent", col = c("pink", "orange"), 
 names.arg = c("not pale", "pale"), names.attr = names(eire),
 criteria = NULL, identify = TRUE, cex.lab = 0.8, pch = c(10, 11),
 xlab = "Are majority people are pale ?", ylab = "Percent", 
 axes = TRUE, lablong = "longitude", lablat = "latitude")
}
\keyword{spatial}
\keyword{univar}
