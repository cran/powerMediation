\name{powerLongFull}
\alias{powerLongFull}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power calculation for longitudinal study with 2 time point
}
\description{
Power calculation for testing if mean changes for 2 groups are the same or not for longitudinal study with 2 time point.
}
\usage{
powerLongFull(delta, 
              sigma1, 
              sigma2, 
              n, 
              rho = 0.5, 
              alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta}{
absolute difference of the mean changes between the two groups: \eqn{\delta=|\mu_1 - \mu_2|} where
\eqn{\mu_1} is the mean change over time \eqn{t} in group 1, 
\eqn{\mu_2} is the mean change over time \eqn{t} in group 2. 
}
  \item{sigma1}{
the standard deviation of baseline values within a treatment group
}
  \item{sigma2}{
the standard deviation of follow-up values within a treatment group
}
  \item{n}{
sample size per group
}
  \item{rho}{
correlation coefficient between baseline and follow-up values within a treatment group.
}
  \item{alpha}{
Type I error rate.
}
}
\details{
The power formula is based on Equation 8.31 on page 336 of Rosner (2006).
\deqn{
power=\Phi\left(-Z_{1-\alpha/2}+\frac{\delta\sqrt{n}}{\sigma_d \sqrt{2}}\right)
}{%
power=\Phi\left(-Z_{1-\alpha/2}+\frac{\delta\sqrt{n}}{\sigma_d \sqrt{2}}\right)
}
where \eqn{\sigma_d = \sigma_1^2+\sigma_2^2-2\rho\sigma_1\sigma_2}, \eqn{\delta=|\mu_1 - \mu_2|},
\eqn{\mu_1} is the mean change over time \eqn{t} in group 1, 
\eqn{\mu_2} is the mean change over time \eqn{t} in group 2, 
\eqn{\sigma_1^2} is the variance of baseline values within a treatment group,
\eqn{\sigma_2^2} is the variance of follow-up values within a treatment group,
\eqn{\rho} is the correlation coefficient between baseline and follow-up values within a treatment group,
and \eqn{Z_u} is the u-th percentile of the standard normal distribution.

We wish to test \eqn{\mu_1 = \mu_2}.

}
\value{
power for testing for difference of mean changes.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ssLong}}, \code{\link{ssLongFull}}, 
\code{\link{powerLong}}.
}
\examples{
    # Example 8.33 on page 336 of Rosner (2006)
    # power=0.80
    powerLongFull(delta=5, sigma1=15, sigma2=15, n=85, rho=0.7, alpha=0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
