\name{dmdirichlet}
\alias{dmdirichlet}
\alias{ddirichlet}
\alias{rdirichlet}
\alias{rmdirichlet}
\alias{modifyDirichlet}
\alias{logIntegrateMixDirichlet}
\alias{logIntegrateDirichlet}

\title{
Functions related to mixture Dirichlet distribution
}
\description{
Functions related to mixture Dirichlet distribution
}
\usage{
dmdirichlet(x, mAlpha, mixtureCoef)
ddirichlet (x, alpha)
rdirichlet (n, alpha) 
rmdirichlet (mAlpha, mixtureCoef) 
modifyDirichlet (prior, y) 
logIntegrateMixDirichlet(y, prior, tau=1) 
logIntegrateDirichlet (y, alpha) 


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector containing a single deviate or matrix containing one random deviate per row.}
  \item{mAlpha}{matrix. Each row is a parameter of Dirichlet}
  \item{alpha}{numeric vector. Parameter for a Dirichlet distribution}
  \item{mixtureCoef}{numeric vector}
  \item{n}{integer}

  \item{prior}{list of two components: alpha and mix.coef}
  \item{y}{numeric vector of counts}
  \item{tau}{numeric}
}
\details{
ddirichlet andn rdirichlet are identically copied from MCMCpack

}
