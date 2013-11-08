\name{powerMediation.VSMc.poisson}
\alias{powerMediation.VSMc.poisson}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for testing mediation effect in poisson regression based on Vittinghoff, Sen and McCulloch's (2009) method
}
\description{
Calculate Power for testing mediation effect in poisson regression based on Vittinghoff, Sen and McCulloch's (2009) method. 
}
\usage{
powerMediation.VSMc.poisson(n, 
                            b2, 
                            sigma.m, 
                            EY, 
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
regression coefficient for the mediator \eqn{m} in the poisson regression
\eqn{\log(E(Y_i))=b0+b1 x_i + b2 m_i}.
}
  \item{sigma.m}{
standard deviation of the mediator.
}
  \item{EY}{
the marginal mean of the outcome.
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
for the poisson regressions:
 \deqn{\log(E(Y_i))=b0+b1 x_i + b2 m_i}{%
\log(E(Y_i))=b0+b1 x_i + b2 m_i
}

Vittinghoff et al. (2009) showed that for the above poisson regression, testing the mediation effect
is equivalent to testing the null hypothesis \eqn{H_0: b_2=0}
versus the alternative hypothesis \eqn{H_a: b_2\neq 0}.

The full model is
 \deqn{\log(E(Y_i))=b_0+b_1 x_i + b_2 m_i }{%
\log(E(Y_i))=b_0+b_1 x_i + b_2 m_i 
 }

The reduced model is
 \deqn{\log(E(Y_i))=b_0+b_1 x_i}{%
\log(E(Y_i))=b_0+b_1 x_i
 }

Vittinghoff et al. (2009) mentioned that if confounders need to be included
in both the full and reduced models, the sample size/power calculation formula
could be accommodated by redefining \code{corr.xm} as the multiple
correlation of the mediator with the confounders as well as the predictor.


}
\value{
  \item{power }{power for testing if \eqn{b_2=0}.}
  \item{delta }{\eqn{b_2\sigma_m\sqrt{(1-\rho_{xm}^2) EY}}}, where
\eqn{\sigma_m} is the standard deviation of the mediator \eqn{m},
\eqn{\rho_{xm}} is the correlation between the predictor \eqn{x}
and the mediator \eqn{m}, and \eqn{EY} is the marginal mean of the 
outcome.
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
\code{\link{minEffect.VSMc.poisson}}, 
\code{\link{ssMediation.VSMc.poisson}} 
}
\examples{
  # example in section 5 (page 546) of Vittinghoff et al. (2009).
  # power = 0.7998578
  powerMediation.VSMc.poisson(n = 1239, b2 = log(1.35), 
    sigma.m = sqrt(0.25 * (1 - 0.25)), EY = 0.5, corr.xm = 0.5,
    alpha = 0.05, verbose = TRUE)
}
\keyword{ test }
