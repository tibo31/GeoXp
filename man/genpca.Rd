\name{genpca}
\alias{genpca}
\title{Generalized Principal Component Analysis (PCA)}
\description{
The function \code{genpca} computes a generalized Principal Component Analysis (PCA).
It calculates the principal components, the coordinates of the variables and 
in these principals components axes and the inertia of these principal components. 
}
\usage{
genpca(data, w = rep(1/nrow(data), length = nrow(data)), 
 m = diag(ncol(data)), center = NULL, reduc = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ matrix \eqn{n \times p}{n x p} }
  \item{w}{vector of size n of weight (by default : \eqn{weight=t(1/n,...,1/n)}{weight=t(1/n,...,1/n)})}
  \item{m}{ matrix \eqn{p \times p}{p x p} (by default : metric=Identity matrix)}
  \item{center}{boolean. if TRUE, centered PCA (by default : center=TRUE)}
  \item{reduc}{boolean. if TRUE, reduced PCA (by default : reduce=TRUE)}
}
\details{
Let \deqn{W=diag(w)}{W=diag(w)}
\deqn{x=data=(x_1',...,x_n')'}{x=data=(x_1',...,x_n')'} with \deqn{x_i=(x_i^1,...,x_i^p)}{x_i=(x_i^1,...,x_i^p)} \cr
Let \deqn{1_n=(1,...,1)' }{1_n=(1,...,1)' } with n rows and : \deqn{1_p=(1,...,1)'}{1_p=(1,...,1)'} with p rows.
Normalization of weight : \deqn{w_i=\frac{w_i}{\sum_iw_i}}{w_i=w_i/sum_iw_i}
Vector of means : \deqn{\bar{x}=(\bar{x^1},...,\bar{x^p})' }{bar(x)=(bar(x^1),...,bar(x^p))' } with:
\deqn{\bar{x^j}=\sum_iw_ix_i^j}{bar(x^j)=sum_iw_ix_i^j} \cr
If center=True, \deqn{x_c=x-1_n\bar{x}'}{x_c=x-bar(x)'}
Standart deviation : \deqn{(\sigma^j)^2=\sum_iw_i(x_i^j)^2-(\bar{x^j})^2}{sigma^j^2=sum_iw_i(x_i^j)^2-(bar(x^j))^2}
\deqn{\Sigma=diag((\sigma^1)^2,...,(\sigma^p)^2)'}{Sigma=diag((sigma^1)^2,...,(sigma^p)^2)'}
If reduc=True : \deqn{x_{cr}=x_c \times \Sigma^{-1/2}}{x_{cr}=x_c Sigma^(-1/2)}
Variance-Covariance matrix: \cr
\deqn{C=x_{cr}'Wx_{cr}}{C=x_cr'Wx_cr}
Cholesky decomposition : \eqn{M=LL'}{M=LL'} where M=m \cr
Let \deqn{C_l=LCL'}{C_l=LCL'}
Let U and D as : \deqn{C_lU=UD}{C_lU=UD} with \eqn{D=diag(\lambda_1,...,\lambda_p)}{D=diag(lambda_1,...,lambda_p)} \cr
Let \deqn{V=L'U}{V=L'U} \cr
Then : \cr
Coordinates of individuals in the principals components basis :
\deqn{CC=x_{cr}V}{CC=x_cr V}
Coordinates of variables in principals components :
\deqn{VC=CVD^{-1/2}}{VC=CVD^(-1/2)}
Inertia :
\deqn{I=D1_p}{I=D1_p}

}


\value{
Returns `inertia' vector of size p with percent of inertia of each component (corresponding to I),
`casecoord' matrix \eqn{n \times p}{n x p} (corresponding to matrix CC),
`varcoord' matrix \eqn{p \times n}{p x n} (corresponding to matrix VC0).}

\references{Thibault Laurent, Anne Ruiz-Gazen, Christine Thomas-Agnan (2012), GeoXp: An R Package for Exploratory Spatial Data Analysis. \emph{Journal of Statistical Software}, 47(2), 1-23. \cr \cr

Caussinus H., Fekri M., Hakam S., Ruiz-Gazen A. (2003) , \emph{A monitoring display of Multivariate Outliers}, Computational Statistics and Data Analysis, vol. 44, 1-2, 237-252.}

\author{Thomas-Agnan C., Aragon Y., Ruiz-Gazen A., Laurent T., Robidou L.}

\keyword{spatial}
\keyword{multivariate}
\seealso{\code{\link{clustermap}},\code{\link{pcamap}}}

