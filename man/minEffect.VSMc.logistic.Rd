\name{minEffect.VSMc.logistic}
\alias{minEffect.VSMc.logistic}
\title{
Minimum detectable slope for mediator in logistic regression based on Vittinghoff, Sen and McCulloch's (2009) method
}
\description{
 Calculate minimal detectable slope for mediator given sample size and power in logistic regression based on Vittinghoff, Sen and McCulloch's (2009) method.
}
\usage{
minEffect.VSMc.logistic(n, power, sigma.m, p, corr.xm, alpha = 0.05, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{power}{
power for testing \eqn{b_2=0} for the logistic regression \eqn{\log(p_i/(1-p_i))=b0+b1 x_i + b2 m_i}.
}
  \item{sigma.m}{
standard deviation of the mediator.
}
  \item{p}{
the marginal prevalence of the outcome.
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
for the logistic regressions:
 \deqn{\log(p_i/(1-p_i))=b_0+b_1 x_i + b_2 m_i }{%
\log(p_i/(1-p_i))=b_0+b_1 x_i + b_2 m_i 
 }

Vittinghoff et al. (2009) showed that for the above logistic regression, testing the mediation effect
is equivalent to testing the null hypothesis \eqn{H_0: b_2=0}
versus the alternative hypothesis \eqn{H_a: b_2\neq 0}, if the
correlation \code{corr.xm} between the primary predictor and mediator is non-zero.

The full model is
 \deqn{\log(p_i/(1-p_i))=b_0+b_1 x_i + b_2 m_i }{%
\log(p_i/(1-p_i))=b_0+b_1 x_i + b_2 m_i 
 }

The reduced model is
 \deqn{\log(p_i/(1-p_i))=b_0+b_1 x_i}{%
\log(p_i/(1-p_i))=b_0+b_1 x_i
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
\code{\link{powerMediation.VSMc.logistic}},
\code{\link{ssMediation.VSMc.logistic}} 
}
\examples{
  # example in section 4 (page 545) of Vittinghoff et al. (2009).
  # minimum effect is log(1.5)= 0.4054651 

  minEffect.VSMc.logistic(n = 255, power = 0.8, sigma.m = 1, 
    p = 0.5, corr.xm = 0.5, alpha = 0.05, verbose = TRUE)

}
\keyword{ test }
