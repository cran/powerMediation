\name{minEffect.VSMc}
\alias{minEffect.VSMc}
\title{
Minimum detectable slope
}
\description{
 Calculate minimal detectable slope given sample size and power for simple linear regression.
}
\usage{
minEffect.VSMc(n, power, sigma.m, sigma.e, corr.xm, alpha = 0.05, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{power}{
power for testing \eqn{b_2=0} for the linear regression \eqn{y_i=b0+b1 x_i + b2 m_i + \epsilon_i, \epsilon_i\sim N(0, \sigma_e^2)}.
}
  \item{sigma.m}{
standard deviation of the mediator.
}
  \item{sigma.e}{
standard deviation of the random error term in the linear regression
\eqn{y_i=b0+b1 x_i + b2 m_i + \epsilon_i, \epsilon_i\sim N(0, \sigma_e^2)}.
}
  \item{corr.xm}{
correlation between the predictor \eqn{x} and the mediator \eqn{m}.
}
  \item{alpha}{
type I error rate.
}
  \item{verbose}{
logical. \code{TRUE} means printing minimum absolute detectable effect; \code{FALSE} means not printing minimum absolute detectable effect.
}
}
\details{
The test is for testing the null hypothesis \eqn{b_2=0}
versus the alternative hypothesis \eqn{b_2\neq 0}
for the linear regressions:
 \deqn{y_i=b_0+b_1 x_i + b_2 m_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{e})}{%
y_i=b_0+b_1 x_i + b_2 m_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{e})
 }

Vittinghoff et al. (2009) showed that for the above linear regression, testing the mediation effect
is equivalent to testing the null hypothesis \eqn{H_0: b_2=0}
versus the alternative hypothesis \eqn{H_a: b_2\neq 0}.
}
\value{
  \item{b2 }{minimum absolute detectable effect.}
  \item{res.uniroot }{results of optimization to find the optimal sample size.}
}
\references{
  Vittinghoff, E. and Sen, S. and McCulloch, C.E..
  Sample size calculations for evaluating mediation.
  \emph{Statistics In Medicine}.  2009;28:541-557.
}
\author{
Weiliang Qiu \email{stwxq@channing.harvard.edu}
}
\note{
The test is a two-sided test. Code for one-sided tests will be added later.
}

\seealso{
\code{\link{powerMediation.VSMc}},
\code{\link{ssMediation.VSMc}} 
}
\examples{
  minEffect.VSMc(n=100, power=0.8, sigma.m=0.1, sigma.e=0.2, corr.xm=0.5, 
    alpha = 0.05, verbose = TRUE)
}
\keyword{ test }
