\name{ss.SLR.rho}
\alias{ss.SLR.rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size for testing slope for simple linear regression based on R2
}
\description{
Calculate sample size for testing slope for simple linear regression based on R2.
}
\usage{
ss.SLR.rho(power, 
           rho2, 
           n.lower = 2.01, 
           n.upper = 1e+30, 
           alpha = 0.05, 
           verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
power.
}
  \item{rho2}{
square of the correlation between the outcome and the predictor.
}
  \item{n.lower}{
lower bound of the sample size.
}
  \item{n.upper}{
upper bound o the sample size.
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
The test is a two-sided test. For one-sided tests, please double the 
significance level. For example, you can set \code{alpha=0.10}
to obtain one-sided test at 5\% significance level.
}

\seealso{
\code{\link{minEffect.SLR}},
\code{\link{power.SLR}},
\code{\link{power.SLR.rho}},
\code{\link{ss.SLR}}.
}

\examples{
  ss.SLR.rho(power=0.8, rho2=0.6, alpha = 0.05, verbose = TRUE)

}
\keyword{ test }

