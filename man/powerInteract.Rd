\name{powerInteract}
\alias{powerInteract}
\title{
Power for detecting interaction effect in 2-way ANOVA
}
\description{
Power for detecting interaction effect in 2-way ANOVA.
}
\usage{
powerInteract(nTotal, a, b, effsize, alpha = 0.05, nTests = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nTotal}{
number of observations in total
}
  \item{a}{
number of levels in factor 1
}
  \item{b}{
number of levels in factor 2
}
  \item{effsize}{
effect size
}
  \item{alpha}{
type I error rate
}
  \item{nTests}{
number of tests. (Family-wise type I error rate will be controlled based on Bonferroni's correction)
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
power
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
Weiliang Qiu \email{stwxq@channing.harvard.edu}
}
\note{
The test is a two-sided test. For one-sided tests, please double the 
significance level. For example, you can set \code{alpha=0.10}
to obtain one-sided test at 5\% significance level.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  powerInteract(nTotal=50, a=2, b=3, effsize=1.5, alpha=0.05, nTests=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ test }
