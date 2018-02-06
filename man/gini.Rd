\name{gini}
\alias{gini}
\title{Calculates a Gini Index }
\description{
The function \code{gini()} calculates the Gini Index associated to the variable \code{var}.
}
\usage{
gini(var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{var}{a vector of numerical values of size n}
}
\details{
Let \eqn{x_k}{x_k}, k=1,...,K be the distinct values taken by \code{var}.
For each site, Gini.r returns two pairs of frequencies. The pair (\code{f},\code{g}) where f represents
\deqn{f_k=\frac{1}{n}\sum_i 1(Var_i=x_k)}{f_k=frac(1)(n)sum_i 1(Var_i=x_k)} (n is the length of \code{var}) and g represents
\deqn{g_k=\frac{x_k}{\bar{x}}f_k}{g_k=frac(x_k)(bar(x))f_k}  The pair (F,G) represents the corresponding cumulative frequencies. \cr
The Gini Index is calculated as :
\deqn{I_G=\frac{1}{2\bar{x}}\sum_{i=1}^K\sum_{j=1}^K|x_i-x_j|}{I_G=frac(1)(2bar(x))sum_(i=1)^K sum_(j=1)^K abs(x_i-x_j)}
}

\note{
This function is used in ginimap.R but it is not an interactive function. }

\value{
(f, F, g, G, gini) where f, F, g and G are vectors of size K and gini a numeric value.
}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\seealso{\code{\link{ginimap}}}
\keyword{spatial}
