\name{nonormmoran}
\alias{nonormmoran}
\title{Detection of spatial autocorrelation}
\description{
The function \code{nonnormmoran} is used to detect spatial autocorrelation in
the residuals \eqn{u}{u} from the least squares model \eqn{Y=\beta \times X+u}{Y=beta X+u}.
It calculates Moran's I statistic of the residuals based on the gaussian asymptotic distribution and give a p-value associated to the test of spatial autocorrelation (gaussian version). }
\usage{
nonormmoran(y, x, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{vector of size n of dependent variable}
  \item{x}{matrix \eqn{n \times p}{n x p} containing explanatory variables}
  \item{W}{spatial weight matrix}
}
\details{
W is supposed standartized :
\deqn{I=(n/s)\frac{u'Wu}{u'u}}{I=(n/s)frac(u'Wu)(u'u)}
\deqn{I \sim N(E(I),var(I))}{I sim N(E(I),var(I))}
let \eqn{M=(I-X(X'X)^{-1}X')}{M=(I-X(X'X)^(-1)X')}
\deqn{E(I)=(n/s)\frac{tr(MW)}{n-k}}{E(I)=(n/s)frac(tr(MW))(n-k)}
\deqn{d=\frac{n-p}{n+p+2}}{d=frac(n-p)(n+p+2)}
\deqn{V(I)=(n/s)^2[tr(MWMW')+tr(MW)^2+(tr(MW))^2]/d-E(I)^2}{V(I)=(n/s)^2[tr(MWMW')+tr(MW)^2+(tr(MW))^2]/d-E(I)^2}
\deqn{Z_I=\frac{I-E(I)}{v(I)^{1/2}}}{Z_I=frac(I-E(I))(v(I)^(1/2))}

When W is row-normalised, \eqn{s=n}{s=n} else \eqn{s}{s} is the number of nonzero links}

\note{
This function is closed to the function \code{lm.morantest} included in the spdep package.
nonnormmoran is less confortable to use because it didn't take into account S4 classes.
}

\value{
a list with \code{nobs} is the number of observations,
\code{nvar}, the number of explanatory variables, \code{morani} is the Moran's I statistic estimate,
\code{imean} is E(I), \code{ivar} is var(I), \code{istat} is the normalized Moran's I statistic
(corresponding to \eqn{Z_I}{Z_i}), and \code{prob} the associated p-value.}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Translated into R from Jim Lessage's \emph{Spatial Econometrics Toolbox}, \url{http://www.spatial-econometrics.com/} }

\keyword{spatial}

\examples{
###
# data baltimore
data(baltimore)

# a spatial weight matrix constructed in the matrix format
W <- makeneighborsw(cbind(baltimore$X,baltimore$Y),method="neighbor",4)

# when W is not row-normalised ...
nonormmoran(baltimore$PRICE,cbind(rep(1,nrow(baltimore)),baltimore[,14:15]),W)
# when W is row_normalised ...
nonormmoran(baltimore$PRICE,cbind(rep(1,nrow(baltimore)),baltimore[,14:15]),normw(W))

# If we compare to the function lm.morantest
baltimore.lm<-lm(PRICE~LOTSZ+SQFT,data=baltimore)

lm.morantest(baltimore.lm, mat2listw(W))
}

\seealso{\code{\link{moranplotmap}},\code{\link{neighbourmap}},\code{\link{makeneighborsw}},\code{\link{normw}} }
