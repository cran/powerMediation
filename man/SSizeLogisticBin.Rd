\name{SSizeLogisticBin}
\alias{SSizeLogisticBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating sample size for simple logistic regression with binary predictor
}
\description{
Calculating sample size for simple logistic regression with binary predictor.
}
\usage{
SSizeLogisticBin(p1, 
                 p2, 
                 B, 
                 alpha = 0.05, 
                 power = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
\eqn{pr(diseased | X = 0)}, i.e. the event rate at \eqn{X = 0} in logistic regression
\eqn{logit(p) = a + b X}, where \eqn{X} is the binary predictor.
}
  \item{p2}{
\eqn{pr(diseased | X = 1)}, the event rate at \eqn{X = 1} in logistic regression
\eqn{logit(p) = a + b X}, where \eqn{X} is the binary predictor.
}
  \item{B}{
\eqn{pr(X = 1)}, i.e. proportion of the sample with \eqn{X = 1}
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
\log(p / (1 - p)) = \beta_0 + \beta_1 X
}{%
\log(p / (1 - p)) = \beta_0 + \beta_1 X
}
where \eqn{p = prob(Y = 1)}, \eqn{X} is the binary predictor, 
\eqn{p_1 = pr(diseased | X = 0)},
\eqn{p_2 = pr(diseased| X = 1)}, \eqn{B = pr(X = 1)}, and 
\eqn{p = (1 - B) p_1 + B p_2}.
The sample size formula we used for testing if 
\eqn{\beta_1 = 0}, is Formula (2) in Hsieh et al. (1998):
\deqn{
n = (Z_{1-\alpha/2}[p(1-p)/B]^{1/2} + Z_{power}[p_1(1-p_1)+p_2(1-p_2)(1-B)/B]^{1/2})^2/[ (p_1-p_2)^2 (1-B) ]
}{%
n=(Z_{1-\alpha/2}[p(1-p)/B]^{1/2} + Z_{power}[p_1(1-p_1)+p_2(1-p_2)(1-B)/B]^{1/2})^2/[ (p_1-p_2)^2 (1-B) ]
}
where \eqn{n} is the required total sample size and \eqn{Z_u} is the \eqn{u}-th
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
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{powerLogisticBin}}
}
\examples{
    ## Example in Table I Design (Balanced design with high event rates) 
    ## of Hsieh et al. (1998 )
    ## the sample size is 1281
    SSizeLogisticBin(p1 = 0.4, p2 = 0.5, B = 0.5, alpha = 0.05, power = 0.95)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
