\name{minEffect.VSMc.cox}
\alias{minEffect.VSMc.cox}
\title{
Minimum detectable slope for mediator in cox regression based on Vittinghoff, Sen and McCulloch's (2009) method
}
\description{
 Calculate minimal detectable slope for mediator given sample size and power in cox regression based on Vittinghoff, Sen and McCulloch's (2009) method.
}
\usage{
minEffect.VSMc.cox(n, 
                   power, 
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
  \item{power}{
power for testing \eqn{b_2=0} for the cox regression \eqn{\log(\lambda)=\log(\lambda_0)+b1 x_i + b2 m_i}, where \eqn{\lambda} is the hazard function and
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
logical. \code{TRUE} means printing minimum absolute detectable effect; \code{FALSE} means not printing minimum absolute detectable effect.
}
}
\details{
The test is for testing the null hypothesis \eqn{b_2=0}
versus the alternative hypothesis \eqn{b_2\neq 0}
for the cox regressions:
 \deqn{\log(\lambda)=\log(\lambda_0)+b_1 x_i + b_2 m_i }{%
\log(\lambda)=\log(\lambda_0)+b_1 x_i + b_2 m_i 
 }

Vittinghoff et al. (2009) showed that for the above cox regression, testing the mediation effect
is equivalent to testing the null hypothesis \eqn{H_0: b_2=0}
versus the alternative hypothesis \eqn{H_a: b_2\neq 0}, if the
correlation \code{corr.xm} between the primary predictor and mediator is non-zero.

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
\code{\link{powerMediation.VSMc.cox}},
\code{\link{ssMediation.VSMc.cox}} 
}
\examples{
  # example in section 6 (page 547) of Vittinghoff et al. (2009).
  # minimum effect is = log(1.5) = 0.4054651

  minEffect.VSMc.cox(n = 1399, power = 0.7999916, 
    sigma.m = sqrt(0.25 * (1 - 0.25)), psi = 0.2, corr.xm = 0.3, 
    alpha = 0.05, verbose = TRUE)

}
\keyword{ test }
