\name{ssMediation.Sobel}
\alias{ssMediation.Sobel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size for testing mediation effectd (Sobel's test)
}
\description{
Calculate sample size for testing mediation effect based on Sobel's test.
}
\usage{
ssMediation.Sobel(power, 
                  theta.1a, 
                  lambda.a, 
                  sigma.x, 
                  sigma.m,
                  sigma.epsilon, 
                  n.lower = 1, 
                  n.upper = 1e+30, 
                  alpha = 0.05, 
                  verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
power of the test.
}
  \item{theta.1a}{
regression coefficient for the predictor in the linear regression linking
the predictor \eqn{x} to the mediator \eqn{m} (\eqn{m_i=\theta_0+\theta_{1a} x_i + e_i, e_i\sim N(0, \sigma^2_e)}).
}
  \item{lambda.a}{
regression coefficient for the mediator in the linear regression linking
the predictor \eqn{x} and the mediator \eqn{m} to the outcome \eqn{y} (\eqn{y_i=\gamma+\lambda_a m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}).
}
  \item{sigma.x}{
standard deviation of the predictor.
}
  \item{sigma.m}{
standard deviation of the mediator.
}
  \item{sigma.epsilon}{
standard deviation of the random error term 
in the linear regression linking
the predictor \eqn{x} and the mediator \eqn{m} to the outcome \eqn{y} (\eqn{y_i=\gamma+\lambda_a m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}).
}
  \item{n.lower}{
lower bound of the sample size.
}
  \item{n.upper}{
upper bound of the sample size.
}
  \item{alpha}{
type I error rate.
}
  \item{verbose}{
logical. \code{TRUE} means printing power; \code{FALSE} means not printing power.
}
}
\details{
The sample size is for testing the null hypothesis \eqn{\theta_1\lambda=0}
versus the alternative hypothesis \eqn{\theta_{1a}\lambda_a\neq 0}
for the linear regressions:
 \deqn{m_i=\theta_0+\theta_{1a} x_i + e_i, e_i\sim N(0, \sigma^2_e)}{%
             m_i=\theta_0+\theta_{1a} x_i + e_i, e_i\sim N(0, \sigma^2_e)}
 \deqn{y_i=\gamma+\lambda_a m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}{%
y_i=\gamma+\lambda_a m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})
  }

Test statistic is based on Sobel's (1982) test:
\deqn{Z=\frac{\hat{\theta}_{1a}\hat{\lambda_a}}{\hat{\sigma}_{\theta_{1a}\lambda_a}}
}{%
Z=\frac{\hat{\theta}_{1a}\hat{\lambda_a}}{\hat{\sigma}_{\theta_{1a}\lambda_a}}
}
where \eqn{\hat{\sigma}_{\theta_{1a}\lambda_a}} is the estimated standard deviation
of the estimate \eqn{\hat{\theta}_{1a}\hat{\lambda_a}} using multivariate
delta method: 
\deqn{\sigma_{\theta_{1a}\lambda_a}=\sqrt{\theta_{1a}^2\sigma_{\lambda_a}^2+\lambda_a^2\sigma_{\theta_{1a}}^2}}{%
\sigma_{\theta_{1a}\lambda_a}=\sqrt{\theta_{1a}^2\sigma_{\lambda_a}^2+\lambda_a^2\sigma_{\theta_{1a}}^2}
}
and
 \eqn{\sigma_{\theta_{1a}}^2=\sigma_e^2/(n\sigma_x^2)} is the 
variance
of the estimate \eqn{\hat{\theta}_{1a}}, and
\eqn{\sigma_{\lambda_a}^2=\sigma_{\epsilon}^2/(n\sigma_m^2(1-\rho_{mx}^2))}
is the variance
of the estimate \eqn{\hat{\lambda_a}}, \eqn{\sigma_m^2} is the variance of the
mediator \eqn{m_i}.

From the linear regression \eqn{m_i=\theta_0+\theta_{1a} x_i+e_i}, we have the
relationship \eqn{\sigma_e^2=\sigma_m^2(1-\rho^2_{mx})}. Hence, we can simply
the variance \eqn{\sigma_{\theta_{1a}, \lambda_a}} to 
\deqn{\sigma_{\theta_{1a}\lambda_a}=\sqrt{\theta_{1a}^2\frac{\sigma_{\epsilon}^2}{n\sigma_m^2(1-\rho_{mx}^2)}+\lambda_a^2\frac{\sigma_{m}^2(1-\rho_{mx}^2)}{n\sigma_x^2}}}{%
\sigma_{\theta_{1a}\lambda_a}=\sqrt{\theta_{1a}^2\frac{\sigma_{\epsilon}^2}{n\sigma_m^2(1-\rho_{mx}^2)}+\lambda_a^2\frac{\sigma_{m}^2(1-\rho_{mx}^2)}{n\sigma_x^2}}}


}
\value{
  \item{n }{sample size.}
  \item{res.uniroot }{results of optimization to find the optimal sample size.}
}
\references{
  Sobel, M. E.
  Asymptotic confidence intervals for indirect effects in structural equation models.
  \emph{Sociological Methodology}. 1982;13:290-312.
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
\code{\link{powerMediation.Sobel}}, 
\code{\link{testMediation.Sobel}} 
}
\examples{
 ssMediation.Sobel(power=0.8, theta.1a=0.1701, lambda.a=0.1998, 
   sigma.x=0.57, sigma.m=0.61, sigma.epsilon=0.2, 
   alpha = 0.05, verbose = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
