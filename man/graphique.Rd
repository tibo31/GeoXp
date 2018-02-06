\name{graphique}
\alias{graphique}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Drawing a graphic}
\description{
The function `graphique' draws a specific graphic. Possibility to draw
one of these graphics : Histogram, Barplot, Boxplot, Polyboxplot, Scatterplot, Densityplot, Neighbourplot, Lorentz, Moran, Angleplot, Variocloud, Cluster and PCA (graphics of individuals and variables). This function is called in most of the functions of GeoXp (this is not an interactive function). 
}

\usage{
graphique(var1, var2, var3, obs, num, graph = "", couleurs = "",
    symbol = 16, labvar = "", nbcol = 10, alpha1, W,
    Xpoly, Ypoly, F, G, opt1 = 1, opt2 = 1, quantiles = 0,
    labmod = "", direct, inertie, label = 0, kernel, obsq, locmoran = FALSE,
    bin = NULL, cex.lab=1, buble=FALSE, cbuble=NULL, legmap=NULL, 
    legends=list(FALSE,FALSE), xlim, ylim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var1}{1st variable of size n}
  \item{var2}{2nd variable of size n}
  \item{var3}{3rd variable of size n (used for variocloudmap)}
  \item{obs}{a vector of boolean of size n with sites selected}
  \item{num}{number of windows which must be activated}
  \item{graph}{name of graphic which must be drawn : Histogram, Barplot, Boxplot, Polyboxplot, Scatterplot, Densityplot1, Densityplot2, Neighbourplot, Lorentz, Moran, Quadrant, Angleplot, Variocloud, Cluster, Acp1, Acp2}
  \item{couleurs}{Possibilty to change colors on the graphic}
  \item{symbol}{choice of representation of selected sites}
  \item{labvar}{name(s) of variable(s) studied}
  \item{nbcol}{number of cells if the graphic choosen is the histogram (10 by default)  }
  \item{alpha1}{regression smoothing paramater}
  \item{W}{Spatial weight matrix}
  \item{Xpoly}{x-coordinates of the vertices of selected polygon}
  \item{Ypoly}{y-coordinates of the vertices of selected polygon}
  \item{F}{Used for Ginimap}
  \item{G}{Used for Ginimap}
  \item{opt1}{Option for adding a curve (regression line or conditionnal quantile depending on the function) }
  \item{opt2}{Option for drawing robust empirical variogram}
  \item{quantiles}{vector which contains the values of \eqn{\alpha}{alpha} for conditionnal quantile}
  \item{labmod}{names of factor if the graphic choosed is a barplot}
  \item{direct}{Used for PCA}
  \item{inertie}{Used for PCA}
  \item{label}{Name of sites}
  \item{kernel}{Name of the kernel choosed in densitymap}
  \item{obsq}{Used for clustermap and barmap}
  \item{locmoran}{Print local moran for each site}
  \item{bin}{The bins chosen to calculate empirical variogram}
  \item{cex.lab}{character size of label}
  \item{buble}{For Local Moran plot only : a boolean with TRUE for drawing bubbles, FALSE otherwise}
  \item{cbuble}{For Local Moran plot only : vector of size n with size of each site depending on variable with which bubbles are constructed }
  \item{legmap}{For Local Moran plot only : a list of (numeric,numeric,numeric,numeric,numeric,numeric,character) with the sizes and the corresponding values of the 3 bubbles represented in the legend and the variable name of variable choosed }
  \item{legends}{For Local Moran plot only : a list of c(boolean,boolean,c(numeric,numeric),c(numeric,numeric)) with TRUE for drawing legends of bubbles (1st argument) or factors (2nd argument) and the coordinates of the location of the upper left corner of the legend box}
  \item{xlim}{the x limits of the plot}
  \item{ylim}{the y limits of the plot} 
 
}
\details{This function is called by any function which draws a graphic. A lot of options
are considered because of the large number of graphics proposed
} 
\value{No values, only drawing of a graphic}
\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y.,  Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}
\keyword{utilities}