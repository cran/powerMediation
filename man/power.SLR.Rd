\name{power.SLR}
\alias{power.SLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for testing slope for simple linear regression
}
\description{
Calculate power for testing slope for simple linear regression.
}
\usage{
power.SLR(n, 
          lambda.a, 
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
  \item{lambda.a}{
regression coefficient in the simple linear regression 
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
logical. \code{TRUE} means printing power; \code{FALSE} means not printing power.
}
}
\details{
The power is for testing the null hypothesis \eqn{\lambda=0}
versus the alternative hypothesis \eqn{\lambda\neq 0}
for the simple linear regressions:
 \deqn{y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{e})}{%
y_i=\gamma+\lambda x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{e})
 }
}
\value{
  \item{power }{power for testing if \eqn{b_2=0}.}
  \item{delta }{\eqn{\lambda\sigma_x\sqrt{n}/\sqrt{\sigma_y^2-(\lambda\sigma_x)^2}}.}
  \item{s}{\eqn{\sqrt{\sigma_y^2-(\lambda\sigma_x)^2}}.}
  \item{t.cr}{\eqn{\Phi^{-1}(1-\alpha/2)}, where \eqn{\Phi} is the cumulative distribution function of the standard normal distribution.}
  \item{rho}{correlation between the predictor \eqn{x} and outcome \eqn{y}
\eqn{=\lambda\sigma_x/\sigma_y}.}
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
The test is a two-sided test. For one-sided tests, please double the 
significance level. For example, you can set \code{alpha=0.10}
to obtain one-sided test at 5\% significance level.
}

\seealso{
\code{\link{minEffect.SLR}},
\code{\link{power.SLR.rho}},
\code{\link{ss.SLR.rho}},
\code{\link{ss.SLR}}.
}
\examples{
  power.SLR(n=100, lambda.a=0.8, sigma.x=0.2, sigma.y=0.5, 
    alpha = 0.05, verbose = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
