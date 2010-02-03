\name{ss.SLR}
\alias{ss.SLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size for testing slope for simple linear regression
}
\description{
Calculate sample size for testing slope for simple linear regression.
}
\usage{
ss.SLR(power, lambda.a, sigma.x, sigma.y, n.lower = 2.01, n.upper = 1e+30, alpha = 0.05, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
power for testing if \eqn{\lambda=0} for the simple linear regression 
\eqn{y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma_{e}^2).}
}
  \item{lambda.a}{
regression coefficient in the simple linear regression 
\eqn{y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma_{e}^2).}
}
  \item{sigma.x}{
standard deviation of the predictor.
}
  \item{sigma.y}{
standard deviation of the outcome.
}
  \item{n.lower}{
lower bound for the sample size.
}
  \item{n.upper}{
upper bound for the sample size.
}
  \item{alpha}{
type I error rate.
}
  \item{verbose}{
logical. \code{TRUE} means printing sample size; \code{FALSE} means not printing sample size.
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
 \item{n }{sample size.}
  \item{res.uniroot }{results of optimization to find the optimal sample size.}
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
\code{\link{minEffect.SLR}},
\code{\link{power.SLR}},
\code{\link{power.SLR.rho}},
\code{\link{ss.SLR.rho}}.
}

\examples{
  ss.SLR(power=0.8, lambda.a=0.8, sigma.x=0.2, sigma.y=0.5, 
    alpha = 0.05, verbose = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
