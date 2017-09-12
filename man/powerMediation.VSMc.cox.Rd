\name{powerMediation.VSMc.cox}
\alias{powerMediation.VSMc.cox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for testing mediation effect in cox regression based on Vittinghoff, Sen and McCulloch's (2009) method
}
\description{
Calculate Power for testing mediation effect in cox regression based on Vittinghoff, Sen and McCulloch's (2009) method. 
}
\usage{
powerMediation.VSMc.cox(n, 
                        b2, 
                        sigma.m, 
                        psi, 
                        corr.xm, 
                        alpha = 0.05, 
                        verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{b2}{
regression coefficient for the mediator \eqn{m} in the cox regression
\eqn{\log(\lambda)=\log(\lambda_0)+b1 x_i + b2 m_i}, where \eqn{\lambda} is the hazard function and
\eqn{\lambda_0} is the baseline hazard function.
}
  \item{sigma.m}{
standard deviation of the mediator.
}
  \item{psi}{
the probability that an observation is uncensored, so that
the number of event \eqn{d= n * psi}, where \eqn{n} is the sample size.
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
for the cox regressions:
 \deqn{\log(\lambda)=\log(\lambda_0)+b_1 x_i + b_2 m_i }{%
\log(\lambda)=\log(\lambda_0)+b_1 x_i + b_2 m_i,
 }
where \eqn{\lambda} is the hazard function and
\eqn{\lambda_0} is the baseline hazard function.

Vittinghoff et al. (2009) showed that for the above cox regression, testing the mediation effect
is equivalent to testing the null hypothesis \eqn{H_0: b_2=0}
versus the alternative hypothesis \eqn{H_a: b_2\neq 0}.

The full model is
 \deqn{\log(\lambda)=\log(\lambda_0)+b_1 x_i + b_2 m_i }{%
\log(\lambda)=\log(\lambda_0)+b_1 x_i + b_2 m_i 
 }

The reduced model is
 \deqn{\log(\lambda)=\log(\lambda_0)+b_1 x_i }{%
\log(\lambda)=\log(\lambda_0)+b_1 x_i 
 }

Vittinghoff et al. (2009) mentioned that if confounders need to be included
in both the full and reduced models, the sample size/power calculation formula
could be accommodated by redefining \code{corr.xm} as the multiple
correlation of the mediator with the confounders as well as the predictor.


}
\value{
  \item{power }{power for testing if \eqn{b_2=0}.}
  \item{delta }{\eqn{b_2\sigma_m\sqrt{(1-\rho_{xm}^2) psi}}}, where
\eqn{\sigma_m} is the standard deviation of the mediator \eqn{m},
\eqn{\rho_{xm}} is the correlation between the predictor \eqn{x}
and the mediator \eqn{m}, and \eqn{psi} is 
the probability that an observation is uncensored, so that
the number of event \eqn{d= n * psi}, where \eqn{n} is the sample size.
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
The test is a two-sided test. For one-sided tests, please double the 
significance level. For example, you can set \code{alpha=0.10}
to obtain one-sided test at 5\% significance level.
}

\seealso{
\code{\link{minEffect.VSMc.cox}}, 
\code{\link{ssMediation.VSMc.cox}} 
}
\examples{
  # example in section 6 (page 547) of Vittinghoff et al. (2009).
  # power = 0.7999916
  powerMediation.VSMc.cox(n = 1399, b2 = log(1.5), 
    sigma.m = sqrt(0.25 * (1 - 0.25)), psi = 0.2, corr.xm = 0.3,
    alpha = 0.05, verbose = TRUE)
}
\keyword{ test }
