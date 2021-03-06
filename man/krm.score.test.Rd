\name{krm.score.test}
\alias{krm.score.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjusted Score Test
}
\description{
Adjusted score test for kernel-based regression models. This function is typically not used directly, but is called within krm.most().
}
\usage{
krm.score.test(formula, data, K, regression.type=c("logistic","linear"), verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object. Model under null.
}
  \item{data}{
a data frame
}
  \item{K}{
a n by n kernel/correlation matrix
}
  \item{regression.type}{
a string
}
  \item{verbose}{
Boolean
}
}
\examples{

dat=sim.liu.2008(n=100, a=0, seed=1)
z=as.matrix(subset(dat, select=c(z.1,z.2,z.3,z.4,z.5)))
rho=1
K=kyotil::getK(z,kernel="rbf",para=rho^-2)
krm.score.test (y~x, dat, K, regression.type="logistic") 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
