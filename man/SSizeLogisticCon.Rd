\name{SSizeLogisticCon}
\alias{SSizeLogisticCon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating sample size for simple logistic regression with continuous predictor
}
\description{
Calculating sample size for simple logistic regression with continuous predictor.
}
\usage{
SSizeLogisticCon(p1, 
                 OR, 
                 alpha = 0.05, 
                 power = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
the event rate at the mean of the continuous predictor \code{X} in logistic regression
\eqn{logit(p) = a + b X},
}
  \item{OR}{
Expected odds ratio. \eqn{\log(OR)} is the change in log odds for the difference between at the mean of \eqn{X} and at one SD above the mean.
}
  \item{alpha}{
Type I error rate.
}
  \item{power}{
power for testing if the odds ratio is equal to one.
}
}
\details{
The logistic regression mode is
\deqn{
\log(p/(1-p)) = \beta_0 + \beta_1 X
}{%
\log(p/(1-p)) = \beta_0 + \beta_1 X
}
where \eqn{p=prob(Y=1)}, \eqn{X} is the continuous predictor, and \eqn{\log(OR)} is the
the change in log odds for the difference between at the mean of \eqn{X} and at one SD above the mean.
The sample size formula we used for testing if \eqn{\beta_1=0} or equivalently
\eqn{OR=1}, is Formula (1) in Hsieh et al. (1998):
\deqn{
n=(Z_{1-\alpha/2} + Z_{power})^2/[ p_1 (1-p_1) [log(OR)]^2 ]
}{%
n=(Z_{1-\alpha/2} + Z_{power})^2/[ p_1 (1-p_1) [log(OR)]^2 ]
}
where \eqn{n} is the required total sample size, \eqn{OR} is the
odds ratio to be tested, \eqn{p_1} is the event rate at the mean
of the predictor \eqn{X}, and \eqn{Z_u} is the \eqn{u}-th
percentile of the standard normal distribution.
}
\value{
total sample size required.
}
\references{
    Hsieh, FY, Bloch, DA, and Larsen, MD.
    A SIMPLE METHOD OF SAMPLE SIZE CALCULATION FOR LINEAR AND LOGISTIC REGRESSION.
    \emph{Statistics in Medicine}. 1998; 17:1623-1634.
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
\code{\link{powerLogisticCon}}
}
\examples{
    ## Example in Table II Design (Balanced design (1)) of Hsieh et al. (1998 )
    ## the sample size is 317
    SSizeLogisticCon(p1 = 0.5, OR = exp(0.405), alpha = 0.05, power = 0.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
