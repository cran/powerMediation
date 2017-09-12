\name{ssLong}
\alias{ssLong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample size calculation for longitudinal study with 2 time point
}
\description{
Sample size calculation for testing if mean changes for 2 groups are the same or not for longitudinal study with 2 time point.
}
\usage{
ssLong(es, 
       rho = 0.5, 
       alpha = 0.05, 
       power = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
effect size of the difference of mean change.
}
  \item{rho}{
correlation coefficient between baseline and follow-up values within a treatment group.
}
  \item{alpha}{
Type I error rate.
}
  \item{power}{
power for testing for difference of mean changes.
}
}
\details{
The sample size formula is based on Equation 8.30 on page 335 of Rosner (2006).
\deqn{
n=\frac{2\sigma_d^2 (Z_{1-\alpha/2} + Z_{power})^2}{\delta^2}
}{%
n=\frac{2\sigma_d^2 (Z_{1-\alpha/2} + Z_{power})^2}{\delta^2}
}
where \eqn{\sigma_d = \sigma_1^2+\sigma_2^2-2\rho\sigma_1\sigma_2}, \eqn{\delta=|\mu_1 - \mu_2|},
\eqn{\mu_1} is the mean change over time \eqn{t} in group 1, 
\eqn{\mu_2} is the mean change over time \eqn{t} in group 2, 
\eqn{\sigma_1^2} is the variance of baseline values within a treatment group,
\eqn{\sigma_2^2} is the variance of follow-up values within a treatment group,
\eqn{\rho} is the correlation coefficient between baseline and follow-up values within a treatment group,
and \eqn{Z_u} is the u-th percentile of the standard normal distribution.

We wish to test \eqn{\mu_1 = \mu_2}.

When \eqn{\sigma_1=\sigma_2=\sigma}, then formula reduces to
\deqn{
n=\frac{4(1-\rho)(Z_{1-\alpha/2}+Z_{\beta})^2}{d^2}
}{%
n=\frac{4(1-\rho)(Z_{1-\alpha/2}+Z_{\beta})^2}{d^2}
}
where \eqn{d=\delta/\sigma}.
}
\value{
required sample size per group
}
\references{
    Rosner, B.
    \emph{Fundamentals of Biostatistics}. Sixth edition. Thomson Brooks/Cole. 2006.
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
\code{\link{ssLongFull}}, \code{\link{powerLong}}, 
\code{\link{powerLongFull}}.
}
\examples{
    # Example 8.33 on page 336 of Rosner (2006)
    # n=85
    ssLong(es=5/15, rho=0.7, alpha=0.05, power=0.8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
