\name{powerMediation.Sobel}
\alias{powerMediation.Sobel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power for testing mediation effect  (Sobel's test)
}
\description{
Calculate power for testing mediation effect based on Sobel's test.
}
\usage{
powerMediation.Sobel(n, 
                     theta.1a, 
                     lambda.a, 
                     sigma.x, 
                     sigma.m,
                     rho2.mx, 
                     sigma.epsilon, 
                     alpha = 0.05, 
                     verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size.
}
  \item{theta.1a}{
regression coefficient for the predictor in the linear regression linking
the predictor \eqn{x} to the mediator \eqn{m} (\eqn{m_i=\theta_0+\theta_1 x_i + e_i, e_i\sim N(0, \sigma^2_e)}).
}
  \item{lambda.a}{
regression coefficient for the mediator in the linear regression linking
the predictor \eqn{x} and the mediator \eqn{m} to the outcome \eqn{y} (\eqn{y_i=\gamma+\lambda m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}).
}
  \item{sigma.x}{
standard deviation of the predictor.
}
  \item{sigma.m}{
standard deviation of the mediator.
}
  \item{rho2.mx}{
square of the correlation between the predictor and the mediator.
}
  \item{sigma.epsilon}{
standard deviation of the random error term 
in the linear regression linking
the predictor \eqn{x} and the mediator \eqn{m} to the outcome \eqn{y} (\eqn{y_i=\gamma+\lambda m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}).
}
  \item{alpha}{
type I error.
}
  \item{verbose}{
logical. \code{TRUE} means printing power; \code{FALSE} means not printing power.
}
}
\details{
The power is for testing the null hypothesis \eqn{\theta_1\lambda=0}
versus the alternative hypothesis \eqn{\theta_{1a}\lambda_a\neq 0}
for the linear regressions:
 \deqn{m_i=\theta_0+\theta_1 x_i + e_i, e_i\sim N(0, \sigma^2_e)}{%
             m_i=\theta_0+\theta_1 x_i + e_i, e_i\sim N(0, \sigma^2_e)}
 \deqn{y_i=\gamma+\lambda m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}{%
y_i=\gamma+\lambda m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})
  }

Test statistic is based on Sobel's (1982) test:
\deqn{Z=\frac{\hat{\theta}_1\hat{\lambda}}{\hat{\sigma}_{\theta_1\lambda}}
}{%
Z=\frac{\hat{\theta}_1\hat{\lambda}}{\hat{\sigma}_{\theta_1\lambda}}
}
where \eqn{\hat{\sigma}_{\theta_1\lambda}} is the estimated standard deviation
of the estimate \eqn{\hat{\theta}_1\hat{\lambda}} using multivariate
delta method: 
\deqn{\sigma_{\theta_1\lambda}=\sqrt{\theta_1^2\sigma_{\lambda}^2+\lambda^2\sigma_{\theta_1}^2}}{%
\sigma_{\theta_1\lambda}=\sqrt{\theta_1^2\sigma_{\lambda}^2+\lambda^2\sigma_{\theta_1}^2}
}
and
 \eqn{\sigma_{\theta_1}^2=\sigma_e^2/(n\sigma_x^2)} is the 
variance
of the estimate \eqn{\hat{\theta}_1}, and
\eqn{\sigma_{\lambda}^2=\sigma_{\epsilon}^2/(n\sigma_m^2(1-\rho_{mx}^2))}
is the variance
of the estimate \eqn{\hat{\lambda}}, \eqn{\sigma_m^2} is the variance of the
mediator \eqn{m_i}.

From the linear regression \eqn{m_i=\theta_0+\theta_1 x_i+e_i}, we have the
relationship \eqn{\sigma_e^2=\sigma_m^2(1-\rho^2_{mx})}. Hence, we can simply
the variance \eqn{\sigma_{\theta_1, \lambda}} to 
\deqn{\sigma_{\theta_1\lambda}=\sqrt{\theta_1^2\frac{\sigma_{\epsilon}^2}{n\sigma_m^2(1-\rho_{mx}^2)}+\lambda^2\frac{\sigma_{m}^2(1-\rho_{mx}^2)}{n\sigma_x^2}}}{%
\sigma_{\theta_1\lambda}=\sqrt{\theta_1^2\frac{\sigma_{\epsilon}^2}{n\sigma_m^2(1-\rho_{mx}^2)}+\lambda^2\frac{\sigma_{m}^2(1-\rho_{mx}^2)}{n\sigma_x^2}}}
}
\value{
  \item{power }{power of the test for the parameter \eqn{\theta_1\lambda}}
  \item{delta }{\eqn{\theta_1\lambda/(sd(\hat{\theta}_1)sd(\hat{\lambda}))}}
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
The test is a two-sided test. Code for one-sided tests will be added later.
}

\seealso{
\code{\link{ssMediation.Sobel}}, 
\code{\link{testMediation.Sobel}} 
}
\examples{
 powerMediation.Sobel(n=248, theta.1a=0.1701, lambda.a=0.1998, 
   sigma.x=0.57, sigma.m=0.61, rho2.mx=0.3, sigma.epsilon=0.2, 
   alpha = 0.05, verbose = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
