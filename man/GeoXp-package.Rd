\name{GeoXp-package}
\alias{GeoXp-package}
\alias{GeoXp}
\docType{package}
\title{
Interactive exploratory spatial data analysis
}
\description{
This package is a tool for researchers in spatial statistics, spatial econometrics, 
geography, ecology etc allowing to link dynamically statistical plots with 
elementary maps. This coupling consists in the fact that the selection of a 
zone on the map results in the automatic highlighting of the corresponding points 
on the statistical graph or reversely, the selection of a portion of the graph 
results in the automatic highlighting of the corresponding points on the map. 
GeoXp includes tools from different areas of spatial statistics 
including geostatistics as well as spatial econometrics and point processes. 
Besides elementary plots like boxplots, histograms or simple scatterplos, GeoXp 
also couples with maps Moran scatterplots, variogram cloud, Lorentz Curves,...
In order to make the most of the multidimensionality of the data, GeoXp includes 
some dimension reduction techniques such as PCA.
}
\details{
\tabular{ll}{
Package: \tab GeoXp\cr
Type: \tab Package\cr
Version: \tab 2.0.0 \cr
Date: \tab 2023-04-05\cr
License: \tab GPL Vesion 2 or later \cr
}

In the version 2.0.0, GeoXp has adopted the \code{sf} object proposed by E. Pebezma
in \code{sf} package. The main advantage of using this structure object is on one hand,
a \code{sf} object contains spatial coordinates and behaves as a \code{data.frame}. \cr
On the map, the coordinates of sites are represented by
using the function \code{st_coordinates} and \code{st_point_on_surface} included in \code{sf} package, which calculates
\code{longitude} (for x-axis) and \code{latitude} (for y-axis), applied on a sf Object.  \cr
In GeoXp, we can find three main groups of functions: \cr
- functions using only one variable: the interest variable
is designed by argument \code{name.var}, a character corresponding to a column of the \code{data.frame}
included in \code{sf.obj}, i.e. the Spatial Class object. It can be a numeric variable (\code{histomap()},
\code{densitymap()}, \code{angleplotmap}...) or a factor variable (or character) (\code{barmap()},...).  \cr
- functions using both several variables: the variables of interest are designed by argument
\code{names.var}, a vector of character corresponding to columns of the \code{data.frame}
included in \code{sf.obj}. It can be two numeric variables (\code{dblehistomap}, \code{dbledensitymap}),
one numeric variable and one factor (\code{histobarmap()}, \code{polyboxplotmap()}),
several numeric variables (\code{plot3dmap}, \code{pcamap()} and \code{clustermap()}). \cr
- functions using both a variable and a spatial weight matrix created as a \code{nb} or \code{listw} object
(see package \code{spdep}). \cr
In the case where \code{sf.obj} is a \code{SpatialPolygonDataFrame}, user will have the opportunity to draw
the polygons of Spatial unit by using the \code{Draw Saptial contours} button in the Tk window. User can also give another sf object as background map with option \code{carte}. \cr
Among options which are common to each function, users have the 
possibility to give a \code{criteria}, vector of boolean of size the number of Spatial units,
with TRUE on specific sites. These sites are then represented by a green croice on the map
by clicking on \code{preselected sites} button on the Tk window. \cr
Moreover, users have the possibility to make bubbles and add some graphs 
(histogram, barplot or scattermap). The potential variables are included
in the \code{data.frame} of the \code{SpatialObject}. Users can choose a proportional
symbol mapping: in function plot, we give value
\eqn{var^{0.5}}{sqrt(var)}. User can choose if a legend has to appear on
the map. He could choose then three values represented by bubbles of 
corresponding sizes. \cr
Finally, users can choose to represent the graphical with different colors
using argument \code{col}. In the case of factors (as function \code{barmap}), users could choose
if a legend with corresponding colors will appear on the map. Users can also modify the representation of
selected sites on map with argument \code{pch}.    \cr
Recent functions \code{barnbmap} and \code{histnbmap} give the opportunity to analyse spatial
weight matrix build using functions included in \code{spdep} package.
}

\author{
Christine Thomas-Agnan, Yves Aragon, Anne Ruiz-Gazen, Thibault Laurent, Laurianne Robidou
Maintainer: Thibault Laurent <thibault.laurent@univ-tlse1.fr>
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Roger S.Bivand, Edzer J.Pebesma, Virgilio Gomez-Rubio (2009),  \emph{Applied Spatial Data Analysis with R}, Springer.
}

\keyword{spatial}

