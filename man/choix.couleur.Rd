\name{choix.couleur}
\alias{choix.couleur }
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of colors before plotting map }

\description{
The function \code{choix.couleur()} is used in most of the GeoXp functions to initialize some parameters before plotting map in
the case of the use of a factor. 
}

\usage{
choix.couleur(graphChoice, listvar = NULL, listnomvar = NULL,
 varChoice1 = NULL, legends, col, pch, spdf = FALSE, num_graph, 
 num_carte)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{graphChoice}{kind of graphic chosen}
   \item{listvar}{list of variables}
   \item{listnomvar}{names of variables}
   \item{varChoice1}{the name of the chosen variable}
   \item{legends}{parameters of plot}
   \item{col}{color of plot}
   \item{pch}{symbols}
   \item{spdf}{a boolean}
   \item{num_graph}{number of the device for the figure}
   \item{num_carte}{number of the device for the map}
}

\details{This function is not an interactive function.}

\value{A list of parameters which permit to use the function carte.}
\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}
\author{Laurent T.}

\keyword{spatial}
\keyword{utilities} 
