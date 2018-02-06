\name{immob}
\alias{immob}
\docType{data}
\title{Some price indices of real estate from biggest cities in France}
\description{
This data frame contains price indices of real estate from biggest cities in France in 2008.
}
\usage{data(immob)}
\format{
  A data frame with 94 observations on the following 10 variables.
  \describe{
    \item{\code{Nom}}{Name of the city}
    \item{\code{Code.INSEE}}{Code Insee of the city}
    \item{\code{Code.region}}{Code of the 'region' area of the city}
    \item{\code{longitude}}{A vector x of size n}
    \item{\code{latitude}}{A vector y of size n}
    \item{\code{prix.vente}}{Average selling price at square metres}
    \item{\code{variation.vente}}{Rates of change in the average selling price for a quarter year}
    \item{\code{prix.location}}{ Average price of rental square metres }
    \item{\code{variation.location}}{Rates of change in the average price of rent for a quarter}
    \item{\code{rentabilite}}{Rentability}
  }
}

\source{
Prepared by T. Laurent.
}

\references{
\url{http://www.fnaim.fr/}
}

\examples{
data(immob)
}
\keyword{datasets}
