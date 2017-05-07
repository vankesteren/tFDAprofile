\name{tFDAprofile}
\alias{tFDAprofile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Function to calculate Q-prop and tFDA-score on a high-dimensional dataset.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function inspects all the \eqn{\frac{p^2-p}{2}}{(p² - p)/2} pairs in an \eqn{n \times p}{n * p} dataset sequentially to determine whether they fall in the fourth quadrant (defined by projt_min and maxt_max). If they do, it calculates their tFDA-score.
}
\usage{
tFDAprofile(X, clss, maxt_max=1.96, projt_min=1.96, display_progress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{(High-dimensional) data matrix with \eqn{n}{n} rows and \eqn{p}{p} columns.}
  \item{cl}{Binary class labels - either a factor or numeric (1/0).}
  \item{maxt_max}{Maximum absolute univariate t-score to be labelled a \code{bivariate} case.}
  \item{projt_min}{Minimum absolute projected t-score to be labelled a \code{bivariate} case.}
  \item{display_progress}{Whether to display a progress bar in the console (not yet available, always true).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{total}{The number of pairs investigated.}
\item{qn}{The number of pairs categorised as \code{bivariate} cases.}
\item{qprop}{The Q-prop value \eqn{\frac{qn}{total}}{qn/total}.}
\item{tfda}{The mean tFDA score of the \code{bivariate} cases.}
%% ...
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
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{tFDAprofile-package}}
}
\examples{
# Test the function on a random 20 * 1000 matrix with random class labels
tFDAprofile(matrix(rnorm(20000),20),rbinom(20,1,0.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Q-prop }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~tFDA }% __ONLY ONE__ keyword per line
