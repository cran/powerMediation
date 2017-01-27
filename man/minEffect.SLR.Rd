\name{minEffect.SLR}
\alias{minEffect.SLR}
\title{
Minimum detectable slope
}
\description{
 Calculate minimal detectable slope given sample size and power for simple linear regression.
}
\usage{
minEffect.SLR(n, 
              power, 
              sigma.x, 
              sigma.y, 
              alpha = 0.05, 
              verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{power}{
power for testing if \eqn{\lambda=0} for the simple linear regression 
\eqn{y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma_{e}^2).}
}
  \item{sigma.x}{
standard deviation of the predictor \eqn{sd(x)}.
}
  \item{sigma.y}{
marginal standard deviation of the outcome \eqn{sd(y)}.
(not the marginal standard deviation \eqn{sd(y|x)})
}
  \item{alpha}{
type I error rate.
}
  \item{verbose}{
logical. \code{TRUE} means printing minimum absolute detectable effect; \code{FALSE} means not printing minimum absolute detectable effect.
}

}
\details{
The test is for testing the null hypothesis \eqn{\lambda=0}
versus the alternative hypothesis \eqn{\lambda\neq 0}
for the simple linear regressions:
 \deqn{y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{e})}{%
y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{e})
 }
}
\value{
 \item{lambda.a }{minimum absolute detectable effect.}
  \item{res.uniroot }{results of optimization to find the optimal minimum absolute detectable effect.}
}
\references{
  Dupont, W.D. and Plummer, W.D..
  Power and Sample Size Calculations for Studies Involving Linear Regression.
  \emph{Controlled Clinical Trials}. 1998;19:589-601.
}
\author{
Weiliang Qiu \email{stwxq@channing.harvard.edu}
}
\note{
The test is a two-sided test. Code for one-sided tests will be added later.
}

\seealso{
\code{\link{power.SLR}},
\code{\link{power.SLR.rho}},
\code{\link{ss.SLR}},
\code{\link{ss.SLR.rho}}.
}

\examples{
  minEffect.SLR(n=100, power=0.8, sigma.x=0.2, sigma.y=0.5, 
    alpha = 0.05, verbose = TRUE)
}
\keyword{ test }
