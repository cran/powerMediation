\name{power.SLR.rho}
\alias{power.SLR.rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for testing slope for simple linear regression
}
\description{
Calculate power for testing slope for simple linear regression.
}
\usage{
power.SLR.rho(n, 
              rho2, 
              alpha = 0.05, 
              verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{rho2}{
square of the correlation between the outcome and the predictor.
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
  \item{delta }{\eqn{\sqrt{n}/\sqrt{1/\rho^2-1}}.}
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
\code{\link{ss.SLR.rho}},
\code{\link{ss.SLR}}.
}
\examples{
  power.SLR.rho(n=100, rho2=0.6, alpha = 0.05, verbose = TRUE)

}
\keyword{ test }
