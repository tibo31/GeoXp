\name{choix.bubble}
\alias{choix.bubble}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Selection of a variable before plotting bubbles}

\description{
The function \code{choix.bubble()} is used in most of the GeoXp functions to initialize some parameters before plotting bubbles.
}

\usage{
choix.bubble(buble, listvar, listnomvar, legends, 
                         num_graph, num_carte)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{buble}{boolean}
   \item{listvar}{list of variables}
   \item{listnomvar}{names of variables}
   \item{legends}{parameters of plot}
   \item{num_graph}{number of the device for the figure}
   \item{num_carte}{number of the device for the map}
}

\details{This function is not an interactive function.}

\value{A list of parameters which permit to use the function carte.}
\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}
\author{Laurent T.}

\keyword{spatial}
\keyword{utilities} 
