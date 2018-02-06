\name{mp.school.ps}
\alias{mp.school.ps}
\docType{data}
\title{Midi-pyrennees school agreggated to the pseudo-canton levels}
\description{
This data frame contains some information about schools agreggated to the pseudo-canton levels in Midi-Pyrenees region.
The pseudo-canton is the spatial unit level.
}
\usage{data(mp.school.ps)}
\format{
  A data frame with 155 observations on the following 7 variables.
  \describe{
    \item{\code{longitude}}{x-coordinate of the pseudo-canton}
    \item{\code{latitude}}{y-coordinate of the pseudo-canton}
    \item{\code{name.canton}}{Name of the pseudo-canton}
    \item{\code{rurality.rate}}{Ratio of the number of rural communes in the pseudo-canton to the number of communes}
    \item{\code{Nb.students.per.class}}{Average number of students per class}
    \item{\code{Cost.per.student}}{Average cast per student} 
    \item{\code{Code}}{The mean age of the teachers in the school}       
      }
}

\source{
Prepared by T. Laurent.
}

\note{
The variables \code{Cost.per.student} and  \code{Nb.students.per.class}
have been permuted because of the confidentiality of this data set.
}

\examples{
data(mp.school.ps)
}
\keyword{datasets}
