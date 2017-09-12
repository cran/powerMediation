\name{powerLong.multiTime}
\alias{powerLong.multiTime}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power calculation for testing if mean changes for 2 groups are the
    same or not for longitudinal study with more than 2 time points
}
\description{
Power calculation for testing if mean changes for 2 groups are the
    same or not for longitudinal study with more than 2 time points.
}
\usage{
powerLong.multiTime(es, m, nn, sx2, rho = 0.5, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
effect size
}
  \item{m}{
number of subjects
}
  \item{nn}{
number of observations per subject
}
  \item{sx2}{
within subject variance
}
  \item{rho}{
within subject correlation
}
  \item{alpha}{
type I error rate
}
}
\details{
We are interested in comparing the slopes of the 2 groups \eqn{A} and \eqn{B}:
\deqn{
\beta_{1A} = \beta_{1B}
}{%
\beta_{1A} = \beta_{1B}
}
where
\deqn{
Y_{ijA}=\beta_{0A}+\beta_{1A} x_{jA} + \epsilon_{ijA}, j=1, \ldots, nn; i=1, \ldots, m
}{%
Y_{ijA}=\beta_{0A}+\beta_{1A} x_{jA} + \epsilon_{ijA}, j=1, \ldots, nn; i=1, \ldots, m
}
and
\deqn{
Y_{ijB}=\beta_{0B}+\beta_{1B} x_{jB} + \epsilon_{ijB}, j=1, \ldots, nn; i=1, \ldots, m
}{%
Y_{ijB}=\beta_{0B}+\beta_{1B} x_{jB} + \epsilon_{ijB}, j=1, \ldots, nn; i=1, \ldots, m
}

The power calculation formula is (Equation on page 30 of Diggle et al. (1994)):
\deqn{
power=\Phi\left[
-z_{1-\alpha} + \sqrt{\frac{m nn s_x^2 es^2}{2(1-\rho)}}
\right]
}{%
power=\Phi\left[
-z_{1-\alpha} + \sqrt{\frac{m nn s_x^2 es^2}{2(1-\rho)}}
\right]
}
where \eqn{es=d/\sigma}, \eqn{d} is the meaninful differnce of interest, 
\eqn{sigma^2} is the variance of the random error, 
\eqn{\rho} is the within-subject correlation, and
\eqn{s_x^2} is the within-subject variance.
}
\value{
power
}
\references{
  Diggle PJ, Liang KY, and Zeger SL (1994). 
  Analysis of Longitundinal Data.  page 30.
  Clarendon Press, Oxford
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
\code{\link{ssLong.multiTime}}
}
\examples{
  # power=0.8
  powerLong.multiTime(es=0.5/10, m=196, nn=3, sx2=4.22, rho = 0.5, alpha = 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
