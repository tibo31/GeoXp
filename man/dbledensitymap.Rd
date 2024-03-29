\name{dbledensitymap}
\alias{dbledensitymap}
\title{Double Kernel density estimates and map}
\description{
The function \code{dbledensitymap} plots two kernel density estimates from 2 variables included
in  \code{names.var} computed with \code{bkde}, and a map with sites of coordinates
\code{sf.obj}. Each site is associated to a value of \code{names.var[1]} and
\code{names.var[2]} and the two windows are interactive. }

\usage{
dbledensitymap(sf.obj, names.var, kernel = "triweight", 
  criteria = NULL, carte = NULL, identify = NULL,
  cex.lab = 0.8, pch = 16, col = c("grey", "lightblue3"), xlab = c("", ""), 
  ylab = "", axes = FALSE, lablong = "", lablat = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sf.obj}{object of class sf}
  \item{names.var}{a vector of character of size 2; attribute names or column numbers in attribute table}
  \item{kernel}{Smoothing kernel (see help(bkde) for list of options)}
  \item{criteria}{a vector of boolean of size the number of spatial units, which permit to represent preselected sites with a cross, using the tcltk window}
  \item{carte}{matrix with 2 columns for drawing spatial polygonal contours : x and y coordinates of the vertices of the polygon}
  \item{identify}{if not NULL, the name of the variable for identifying observations on the map}
  \item{cex.lab}{character size of label}
  \item{pch}{16 by default, symbol for selected points}
  \item{col}{c("grey","lightblue3") by default, color of the two density curves}
  \item{xlab}{a list of title for the two x-axis graphics}
  \item{ylab}{a list of title for the two y-axis graphics}
  \item{axes}{a boolean with TRUE for drawing axes on the map}
  \item{lablong}{name of the x-axis that will be printed on the map}
  \item{lablat}{name of the y-axis that will be printed on the map}
}

\details{
The user can choose an interval on the density curve by mouse clicking on the lower and upper boundaries of the interval or by giving directly these values.
The selected sites are then represented on the map in red.
A selection by `points' or `polygon' on the map results in the drawing of the density of the corresponding
sub-distribution on the density plot. Finally, the user can modify the bandwith parameter with a cursor in
the tk window (parameter \eqn{\alpha}{alpha}). \eqn{\alpha}{alpha} is the smoothing parameter for the kernel
smooth : it represents the mean percentage of sample points involved in the local averaging
(example : \eqn{\alpha=20}{alpha=20} means that on average, \eqn{n \times 0.2}{n x 0.2} points are in any
interval of length 2h where h is the usual bandwidth).}


\value{
In the case where user click on \code{save results} button,
a vector of integer is created as a global variable in \code{last.select} object.
It corresponds to the number of spatial units selected just before leaving the Tk window.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer. \cr \cr

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. New York: Springer. \cr \cr

Wand M.P. et Jones M.C. (1995), \emph{Kernel Smoothing}, Chapman \& Hall.
}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{histomap}}, \code{\link{histobarmap}}, \code{\link{scattermap}}, \code{\link{densitymap}}}

\examples{

#########
# data auckland
require(sf)
auckland <- st_read(system.file("shapes/auckland.shp", package="spData")[1])

dbledensitymap(auckland, c("M77_85", "Und5_81"), xlab = c("M77_85", "Und5_81"),
  criteria = (auckland$M77_85 > mean(auckland$M77_85)))

\dontrun{ 
# data eire
require(sf)
eire <- st_read(system.file("shapes/eire.shp", package="spData")[1])

dbledensitymap(eire, c("A", "towns"), kernel = "normal",
  xlab = c("Individuals rate of blood type A", 
  "Surface urbaine"), identify = TRUE)
}
}

\keyword{spatial}
\keyword{smooth} 