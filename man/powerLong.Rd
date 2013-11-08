\name{powerLong}
\alias{powerLong}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power calculation for longitudinal study with 2 time point
}
\description{
Power calculation for testing if mean changes for 2 groups are the same or not for longitudinal study with 2 time point.
}
\usage{
powerLong(es, 
          n, 
          rho = 0.5, 
          alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
effect size of the difference of mean change.
}
  \item{n}{
sample size per group.
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

When \eqn{\sigma_1=\sigma_2=\sigma}, then formula reduces to
\deqn{
power=\Phi\left(-Z_{1-\alpha/2} + \frac{|d|\sqrt{n}}{2\sqrt{1-\rho}}\right)
}{%
power=\Phi\left(-Z_{1-\alpha/2} + \frac{|d|\sqrt{n}}{2\sqrt{1-\rho}}\right)
}
where \eqn{d=\delta/\sigma}.

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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ssLong}}, \code{\link{ssLongFull}}, 
\code{\link{powerLongFull}}.
}
\examples{
    # Example 8.34 on page 336 of Rosner (2006)
    # power=0.75
    powerLong(es=5/15, n=75, rho=0.7, alpha=0.05)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
