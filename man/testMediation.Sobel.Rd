\name{testMediation.Sobel}
\alias{testMediation.Sobel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
P-value and confidence interval for testing mediation effect  (Sobel's test)
}
\description{
Calculate p-value and confidence interval for testing mediation effect based on Sobel's test.
}
\usage{
testMediation.Sobel(theta.1.hat, 
                    lambda.hat, 
                    sigma.theta1, 
                    sigma.lambda, 
                    alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta.1.hat}{
estimated regression coefficient for the predictor in the linear regression linking
the predictor \eqn{x} to the mediator \eqn{m} (\eqn{m_i=\theta_0+\theta_1 x_i + e_i, e_i\sim N(0, \sigma^2_e)}).
}
  \item{lambda.hat}{
estimated regression coefficient for the mediator in the linear regression linking
the predictor \eqn{x} and the mediator \eqn{m} to the outcome \eqn{y} (\eqn{y_i=\gamma+\lambda m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}).
}
  \item{sigma.theta1}{
standard deviation of \eqn{\hat{\theta}_1}
in the linear regression linking
the predictor \eqn{x} to the mediator \eqn{m} (\eqn{m_i=\theta_0+\theta_1 x_i + e_i, e_i\sim N(0, \sigma^2_e)}).
}
  \item{sigma.lambda}{
standard deviation of \eqn{\hat{\lambda}}
in the linear regression linking
the predictor \eqn{x} and the mediator \eqn{m} to the outcome \eqn{y} (\eqn{y_i=\gamma+\lambda m_i+ \lambda_2 x_i + \epsilon_i, \epsilon_i\sim N(0, \sigma^2_{\epsilon})}).
}
  \item{alpha}{
significance level of a test.
}
}
\details{
The test is for testing the null hypothesis \eqn{\theta_1\lambda=0}
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
 \eqn{\hat{\sigma}_{\theta_1}} is the estimated standard deviation
of the estimate \eqn{\hat{\theta}_1}, and
\eqn{\hat{\sigma}_{\lambda}} is the estimated standard deviation
of the estimate \eqn{\hat{\lambda}}.
}
\value{
  \item{pval }{p-value for testing the null hypothesis \eqn{\theta_1\lambda=0}
versus the alternative hypothesis \eqn{\theta_{1a}\lambda_a\neq 0}.}
  \item{CI.low }{Lower bound of the \eqn{100 (1-\alpha)\%} confidence interval for the parameter \eqn{\theta_1\lambda}.}
  \item{CI.upp }{Upper bound of the \eqn{100 (1-\alpha)\%} confidence interval for the parameter \eqn{\theta_1\lambda}.
}
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
\code{\link{powerMediation.Sobel}}, 
\code{\link{ssMediation.Sobel}} 
}
\examples{
  testMediation.Sobel(theta.1.hat=0.1701, lambda.hat=0.1998, 
    sigma.theta1=0.01, sigma.lambda=0.02, alpha=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
