\name{tFDAprofile-package}
\alias{tFDAprofile-package}
\alias{tFDAprofile}
\docType{package}
\title{
\packageTitle{tFDAprofile}
}
\description{
\packageDescription{tFDAprofile}
}
\author{
\packageAuthor{tFDAprofile}

Maintainer: \packageMaintainer{tFDAprofile}
}
\references{

}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
\code{\link[tFDAprofile:tFDAprofile]{tFDAprofile}}

}
\examples{
# Profile a random 20 * 1000 matrix with random class labels for tFDA and Q-prop.
tFDAprofile(matrix(rnorm(20000),20),rbinom(20,1,0.5))
}
