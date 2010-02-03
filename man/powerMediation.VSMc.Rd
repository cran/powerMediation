\name{powerMediation.VSMc}
\alias{powerMediation.VSMc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for testing mediation effect based on Vittinghoff and McCulloch's (2009) method
}
\description{
Calculate Power for testing mediation effect based on Vittinghoff and McCulloch's (2009) method. 
}
\usage{
powerMediation.VSMc(n, b2, sigma.m, sigma.e, corr.xm, alpha = 0.05, 
  verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{b2}{
regression coefficient for the mediator \eqn{m} in the linear regression
\eqn{y_i=b0+b1 x_i + b2 m_i + \epsilon_i, \epsilon_i\sim N(0, \sigma_e^2)}.
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
logical. \code{TRUE} means printing power; \code{FALSE} means not printing power.
}
}
\details{
The power is for testing the null hypothesis \eqn{b_2=0}
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
  \item{power }{power for testing if \eqn{b_2=0}.}
  \item{delta }{\eqn{b_2\sigma_m\sqrt{1-\rho_{xm}}/\sigma_e}, where
\eqn{\sigma_m} is the standard deviation of the mediator \eqn{m},
\eqn{\rho_{xm}} is the correlation between the predictor \eqn{x}
and the mediator \eqn{m}, and \eqn{\sigma_e} is the standard deviation
of the random error term in the linear regression.}
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
\code{\link{minEffect.VSMc}}, 
\code{\link{ssMediation.VSMc}} 
}
\examples{
  powerMediation.VSMc(n=100, b2=0.8, sigma.m=0.1, sigma.e=0.2, corr.xm=0.5, 
    alpha = 0.05, verbose = TRUE)
}
\keyword{ test }
