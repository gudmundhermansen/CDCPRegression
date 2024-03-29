% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cdcp.regression.data.R
\name{cdcp.regression.data}
\alias{cdcp.regression.data}
\title{Pre-process Data for the CDCPRegression Package}
\usage{
cdcp.regression.data(
  y,
  X,
  Z = NULL,
  D = NULL,
  index,
  group,
  index_dummy = FALSE,
  group_dummy = FALSE,
  lag = FALSE
)
}
\arguments{
\item{y}{the response variable, a \eqn{(n*m)} numeric vector, where \eqn{m} is the number of individuals and \eqn{n} is the number of observations per individual.}

\item{X}{a numeric covariate matrix of size \eqn{(n*m) \times p}. The matrix X repent the part of the model specification we believe is affected by a change point.}

\item{Z}{protected covariates, i.e. covariates we believe are constant and are not affected by the change point.}

\item{D}{protected categorical covariates.}

\item{index}{a \eqn{(n*m)}-vector of indexes, e.g. time or year, such that index = 1 represent the first observation for each individual in the panel.}

\item{group}{a vector of length \eqn{(n*m)} representing a potential group structure in the data. If all observations belong to the same group, this should be a \eqn{(n*m)}-vector with the value 1.}

\item{index_dummy}{set \code{index_dummy = TRUE} to include a categorical variable for each index (e.g. year dummy) in the model.}

\item{group_dummy}{set \code{group_dummy = TRUE} to include a categorical variable for each group.}

\item{lag}{set \code{lag = TRUE} to add the lagged response \eqn{y_{i - 1}} to the model.}
}
\value{
A list of data used as input to the various functions in the CDCPRegression Package.
}
\description{
This function pre-process and creates a structured set of data to be used with the various functions in the CDCPRegression
package. The model must be within the class \eqn{y_{i, j} = x_{i, j}^{\rm t} \beta_{\rm L} + x_{i, j}^{\rm t} (\beta_{\rm R} - \beta_{\rm L}) I(t_i < \delta) + \gamma Z_{i ,j} + \alpha D_{i, j} + y_{i - 1, j} + \epsilon_{i, j}}; see Hermansen (2021) for more details.
}
\examples{
# Example 1: Simple Illustration with one Individual 

index <- 1:n
group <- rep(1, n)
 
X <- matrix(1, nrow = n, ncol = 1)
y <- 1.5*(index <= 50) + rnorm(n)

data <- cdcp.regression.data(y = y, X = X, index = index, group = group)


# Example 2: Simple Illustration with m Individuals 

m <- 2
n <- 100 

index <- rep(1:n, times = m)
group <- rep(1:m, each = n)

X <- matrix(1, nrow = n*m)
y <- 1.5*(index <= 50) + rnorm(n = m*n)

data <- cdcp.regression.data(y = y, X = X, index = index, group = group)
}
\references{
Hermansen, G., Knutsen, Carl Henrik & Nygaard, Haavard Mokleiv. (2021). Characterizing and assessing temporal heterogeneity: Introducing a change point framework, with applications on the study of democratization, Political Analysis, 29, 485-504

Cunen, C., Hermansen, G., & Hjort, N. L. (2018). Confidence distributions for change-points and regime shifts. Journal of Statistical Planning and Inference, 195, 14-34.
}
