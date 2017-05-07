\name{tFDA}
\alias{tFDA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to calculate tFDA-score for a single variable pair.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function calculates the tFDA-score for a variable pair (\eqn{n \times 2}{n * 2} matrix) associated with binary class labels.
}
\usage{
tFDA(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Matrix with \eqn{n}{n} rows and 2 columns.}
  \item{cl}{Class labels as numeric vector containing 0 and 1.}
}

\value{
Numeric value representing tFDA-score. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Erik-Jan van Kesteren
}
\note{
%%  ~~further notes~~
This function is used inside the main \code{\link{tFDAprofile}} function.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{tFDAprofile}}, \code{\link{tFDAprofile-package}}
}
\examples{
# Calculate the tFDA score of a random pair of 20 observations
tFDAprofile::tFDA(matrix(rnorm(40),20),rbinom(20,1,0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~tFDA }% use one of  RShowDoc("KEYWORDS")
