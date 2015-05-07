\name{powerPoisson}
\alias{powerPoisson}
\title{
Power calculation for simple Poisson regression
}
\description{
Power calculation for simple Poisson regression.
Assume the predictor is normally distributed.
}
\usage{
powerPoisson(
    beta0, 
    beta1, 
    mu.x1, 
    sigma2.x1, 
    mu.T = 1, 
    phi = 1, 
    alpha = 0.05, 
    N = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta0}{
intercept
}
  \item{beta1}{
slope
}
  \item{mu.x1}{
mean of the predictor
}
  \item{sigma2.x1}{
variance of the predictor
}
  \item{mu.T}{
mean exposure time
}
  \item{phi}{
a measure of over-dispersion
}
  \item{alpha}{
type I error rate
}
  \item{N}{
sample size
}
}
\details{
The simple Poisson regression has the following form:
\deqn{
Pr(Y_i = y_i | mu_i, t_i) = \exp(-\mu_i t_i) (\mu_i t_i)^{y_i}/ (y_i!)
}{%
Pr(Y_i = y_i | mu_i, t_i) = \exp(-\mu_i t_i) (\mu_i t_i)^{y_i}/ (y_i!)
}
where
\deqn{ 
\mu_i=\exp(\beta_0+\beta_1 x_{1i})
}{%
\mu_i=\exp(\beta_0+\beta_1 x_{1i})
}
We are interested in testing the null hypothesis \eqn{\beta_1=0}
versus the alternative hypothesis \eqn{\beta_1 = \theta_1}.
Assume \eqn{x_{1}}  is normally distributed with mean
\eqn{\mu_{x_1}} and variance \eqn{\sigma^2_{x_1}}.
The sample size calculation formula derived by Signorini (1991) is
\deqn{
N=\phi\frac{\left[z_{1-\alpha/2}\sqrt{V\left(b_1 | \beta_1=0\right)}
+z_{power}\sqrt{V\left(b_1 | \beta_1=\theta_1\right)}\right]^2}
{\mu_T \exp(\beta_0) \theta_1^2}
}{%
N=\phi\frac{\left[z_{1-\alpha/2}\sqrt{V\left(b_1 | \beta_1=0\right)}
+z_{power}\sqrt{V\left(b_1 | \beta_1=\theta_1\right)}\right]^2}
{\mu_T \exp(\beta_0) \theta_1^2}
}
where \eqn{\phi} is the over-dispersion parameter (\eqn{=var(y_i)/mean(y_i)}), 
\eqn{\alpha}$ is the type I error rate, 
\eqn{b_1} is the estimate of the slope \eqn{\beta_1},
\eqn{\beta_0} is the intercept,
\eqn{\mu_T} is the mean exposure time, 
\eqn{z_{a}} is the \eqn{100*a}-th lower percentile of 
the standard normal distribution, and \eqn{V\left(b_1|\beta_1=\eta\right)}
is the variance of the estimate \eqn{b_1} given the true slope
\eqn{\beta_1=\eta}. 

The variances are
\deqn{
V\left(b_1 | \beta_1 = 0\right)=\frac{1}{\sigma^2_{x_1}}
}{%
V\left(b_1 | \beta_1 = 0\right)=\frac{1}{\sigma^2_{x_1}}
}
and
\deqn{
V\left(b_1 | \beta_1 = \theta_1\right)=\frac{1}{\sigma^2_{x_1}}
\exp\left[-\left(\theta_1 \mu_{x_1} + \theta_1^2\sigma^2_{x_1}/2\right)\right]
}{%
V\left(b_1 | \beta_1 = \theta_1\right)=\frac{1}{\sigma^2_{x_1}}
\exp\left[-\left(\theta_1 \mu_{x_1} + \theta_1^2\sigma^2_{x_1}/2\right)\right]
}

}
\value{
power
}
\references{
    Signorini D.F. (1991). 
    Sample size for Poisson regression. 
    Biometrika. Vol.78. no.2, pp. 446-50
}
\author{
Weiliang Qiu <stwxq@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
See Also as \code{\link{sizePoisson}}
}
\examples{
# power = 0.8090542
print(powerPoisson(
    beta0 = 0.1, 
    beta1 = 0.5, 
    mu.x1 = 0, 
    sigma2.x1 = 1, 
    mu.T = 1, 
    phi = 1, 
    alpha = 0.05, 
    N = 28))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }
