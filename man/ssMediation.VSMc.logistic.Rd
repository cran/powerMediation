\name{ssMediation.VSMc.logistic}
\alias{ssMediation.VSMc.logistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size for testing mediation effect in logistic regression based on Vittinghoff, Sen and McCulloch's (2009) method
}
\description{
Calculate sample size for testing mediation effect in logistic regression based on Vittinghoff, Sen and McCulloch's (2009) method. 
}
\usage{
ssMediation.VSMc.logistic(power, 
                          b2, 
                          sigma.m, 
                          p, 
                          corr.xm, 
                          n.lower = 1, 
                          n.upper = 1e+30, 
                          alpha = 0.05, 
                          verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
power for testing \eqn{b_2=0} for the logistic regression \eqn{\log(p_i/(1-p_i))=b0+b1 x_i + b2 m_i}.
}
  \item{b2}{
regression coefficient for the mediator \eqn{m} in the logistic regression
\eqn{\log(p_i/(1-p_i))=b0+b1 x_i + b2 m_i}.
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
for the logistic regressions:
 \deqn{\log(p_i/(1-p_i))=b_0+b_1 x_i + b_2 m_i}{%
\log(p_i/(1-p_i))=b_0+b_1 x_i + b_2 m_i
 }

Vittinghoff et al. (2009) showed that for the above logistic regression, testing the mediation effect
is equivalent to testing the null hypothesis \eqn{H_0: b_2=0}
versus the alternative hypothesis \eqn{H_a: b_2\neq 0}.

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
The test is a two-sided test. For one-sided tests, please double the 
significance level. For example, you can set \code{alpha=0.10}
to obtain one-sided test at 5\% significance level.
}

\seealso{
\code{\link{minEffect.VSMc.logistic}}, 
\code{\link{powerMediation.VSMc.logistic}} 
}
\examples{
  # example in section 4 (page 545) of Vittinghoff et al. (2009).
  # n=255

  ssMediation.VSMc.logistic(power = 0.80, b2 = log(1.5), sigma.m = 1, p = 0.5, 
    corr.xm = 0.5, alpha = 0.05, verbose = TRUE)

}
\keyword{ test}
