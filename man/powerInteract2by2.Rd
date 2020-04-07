\name{powerInteract2by2}
\alias{powerInteract2by2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power Calculation for Interaction Effect in 2x2 Two-Way ANOVA Given Effect Sizes
}
\description{
Power calculation for interaction effect in 2x2 two-way ANOVA given effect sizes.
}
\usage{
powerInteract2by2(n, tauBetaSigma, alpha = 0.05, nTests = 1, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
integer. Number of subjects per group.
}
  \item{tauBetaSigma}{
Effect sizes \eqn{\left(\tau\beta\right)_{ij}/\sigma, i=1, \ldots, a, j=1,\ldots, b}, where \eqn{a=b=2} and \eqn{\sigma} is
the standard deviation of random error. Rows are for factor 1 and columns are for factor 2.  Note that \eqn{\sum_{i=1}^a \left(\tau\beta\right)_{ij} = \sum_{j=1}^b \left(\tau\beta\right)_{ij}=0}. We can get
\eqn{\left(\tau\beta\right)_{11}=\theta},
\eqn{\left(\tau\beta\right)_{12}=-\theta},
\eqn{\left(\tau\beta\right)_{21}=-\theta},
\eqn{\left(\tau\beta\right)_{22}=\theta}. So
\code{tauBetaSigma=}\eqn{\theta/\sigma}
}
  \item{alpha}{
family-wise type I error rate.
}
  \item{nTests}{
integer. For high-throughput omics study, 
we perform two-way ANOVA for each of 'nTests' probes.
We use Bonferroni correction to control for family-wise type I error rate.
That is, for each probe, type I error rate would be \code{alpha/nTests}.
}
  \item{verbose}{
logical. Indicating if intermediate results should be printed out.
}
}
\details{
We assume the following model:
\deqn{y_{ijk}=\mu+\tau_i + \beta_j + \left(\tau\beta\right)_{ij} + \epsilon_{ijk},}
where \eqn{i=1,\ldots, a, j=1,\ldots, b, k=1, \ldots, n}, \eqn{\sum_{i=1}^{a}\tau_i = 0}, 
\eqn{\sum_{j=1}^{b}\beta_j = 0}, \eqn{\sum_{i=1}^{a} \left(\tau\beta\right)_{ij} = 0},
\eqn{\sum_{j=1}^{b} \left(\tau\beta\right)_{ij} = 0},
and \eqn{\epsilon_{ijk}\stackrel{i.i.d}{\sim} N\left(0, \sigma^2\right)}. 

The group means are
\deqn{\mu_{ij} = \mu+\tau_i + \beta_j + \left(\tau\beta\right)_{ij}, i=1 \ldots, a,
j=1,\ldots, b.}
Note that \eqn{\mu = \sum_{i=1}^{a}\sum_{j=1}^b \mu_{ij} / (ab)},
\eqn{\tau_i = \sum_{j=1}^b \mu_{ij}/b - \mu}, and 
\eqn{\beta_j = \sum_{i=1}^a \mu_{ij}/a - \mu}.

The null hypothesis \eqn{H_0}: all \eqn{\left(\tau\beta\right)_{ij}, i=1, \ldots, a, j=1,\ldots, b} are equal to zero.
The alternative hypothesis \eqn{H_a}: at least one \eqn{\left(\tau\beta\right)_{ij}} is different from zero.

The F test statistic is \deqn{F=MS_{AB}/MS_{E}\stackrel{H_a}{\sim} F_{(a-1)(b-1), ab(n - 1), ncp},}
where \code{ncp} is the 
non-centrality parameter of the F test statistic:
\deqn{ncp=n\sum_{i=1}^{a}\sum_{j=1}^{b}\left[\frac{\left(\tau\beta\right)_{ij}}{\sigma}\right]^2.}

For the scenario \eqn{a=b=2}, we have 
\eqn{\left(\tau\beta\right)_{11}=\theta},
\eqn{\left(\tau\beta\right)_{12}=-\theta},
\eqn{\left(\tau\beta\right)_{21}=-\theta},
\eqn{\left(\tau\beta\right)_{22}=\theta}. 
Hence,
the non-centrality parameter can be simplified to
\deqn{ncp=4n\left(\frac{\theta}{\sigma}\right)^2.}

The power for testing the null hypothesis \eqn{H_0} versus
the alternative hypothesis \eqn{H_a} is
\deqn{power=Pr\left(F > F_0 | H_a\right),} 
where the rejection region boundary \eqn{F_0} satisfies:
\deqn{Pr\left(F > F_0 | H_0\right) = \alpha/nTests.}
}
\value{
A list with 5 elements:
  \item{power }{the power of the two-way ANOVA test}
  \item{df1 }{the first degree of freedom of the F test statistic (\code{df1=(a-1)(b-1)})}
  \item{df2 }{the second degree of freedom of the F test statistic (\code{df1=a*b(n-1)})}
  \item{F0 }{the rejection region boundary}
  \item{ncp }{the non-centrality parameter}

}
\references{
  Chow SC, Shao J, and Wang H. 
  Sample size calculations in clinical research. 2nd edition.
  Chapman & Hall/CRC. 2008
  
  Montgomery DC.
  Design and Analysis of Experiments. 8th edition.
  John Wiley & Sons. Inc.
}
\author{
Weiliang Qiu \email{weiliang.qiu@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
n = 25
tauBetaSigma = 0.3

# power = 0.8437275
res2 = powerInteract2by2(n = n, tauBetaSigma = tauBetaSigma, 
    alpha = 0.05, nTests = 1, verbose = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

