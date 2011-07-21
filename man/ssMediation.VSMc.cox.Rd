\name{ssMediation.VSMc.cox}
\alias{ssMediation.VSMc.cox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size for testing mediation effect in cox regression based on Vittinghoff, Sen and McCulloch's (2009) method
}
\description{
Calculate sample size for testing mediation effect in cox regression based on Vittinghoff, Sen and McCulloch's (2009) method. 
}
\usage{
ssMediation.VSMc.cox(power, b2, sigma.m, psi, corr.xm, 
  n.lower=1, n.upper=1e+30, alpha = 0.05, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
power for testing \eqn{b_2=0} for the cox regression \eqn{\log(\lambda)=\log(\lambda_0)+b1 x_i + b2 m_i}, where \eqn{\lambda} is the hazard function and
\eqn{\lambda_0} is the baseline hazard function.
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
The test is for testing the null hypothesis \eqn{b_2=0}
versus the alternative hypothesis \eqn{b_2\neq 0}
for the cox regressions:
 \deqn{\log(\lambda)=\log(\lambda_0)+b1 x_i + b2 m_i}{%
\log(\lambda)=\log(\lambda_0)+b1 x_i + b2 m_i
 }

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
  \item{n }{sample size.}
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
\code{\link{minEffect.VSMc.cox}}, 
\code{\link{powerMediation.VSMc.cox}} 
}
\examples{
  # example in section 6 (page 547) of Vittinghoff et al. (2009).
  # n = 1399
  ssMediation.VSMc.cox(power = 0.7999916, b2 = log(1.5), 
    sigma.m = sqrt(0.25 * (1 - 0.25)), psi = 0.2, corr.xm = 0.3,
    alpha = 0.05, verbose = TRUE)

}
\keyword{ test}
