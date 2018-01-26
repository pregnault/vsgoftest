\name{vsgoftest-package}
\alias{vsgoftest-package}
\alias{vsgoftest}
\docType{package}
\title{
\packageTitle{vsgoftest}
}
\description{
\packageDescription{vsgoftest}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{vsgoftest}
\packageIndices{vsgoftest}
}
\author{
\packageAuthor{vsgoftest}

Maintainer: \packageMaintainer{vsgoftest}
}
\references{

Vasicek, O., A test for normality based on sample entropy, \emph{Journal of the Royal Statistical Society,} \bold{38(1)}, 54-59 (1976).

Song, K. S., Goodness-of-fit tests based on Kullback-Leibler discrimination information, \emph{Information Theory, IEEE Transactions on,} \bold{48(5)}, 1103-1117 (2002).

Girardin, V., Lequesne, J. Entropy-based goodness-of-fit tests - a unifying framework. Application to DNA replication. \emph{Accepted in Communications in Statistics: Theory and Methods} (2017). 

Lequesne, J., Regnault, P. Package vsgoftest for R: goodness-of-fit tests based on Kullback-Leibler divergence. \emph{arxiv:} (2018).

}
\keyword{ package }

\examples{
set.seed(1)
samp <- rnorm(50, mean = 2, s = 3)

##Estimating entropy
entropy.estimate(x = samp, window = 8)
log(2*pi*exp(1))/2 #true value of entropy of normal distribution

##Testing normality
vs.test(x = samp, densfun = 'dnorm', param = c(2,3), B = 5000) #Simple null hypothesis
vs.test(x = samp, densfun='dnorm') #Composite null hypothesis

}
