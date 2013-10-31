\name{hmmMargLlik}
\alias{hmmMargLlik}
\alias{readPriorFromFile}
\title{
Functions related to profile HMM}
\description{
Functions related to profile HMM}
\usage{
hmmMargLlik(dat, aaPrior, tau)
readPriorFromFile(priorFileName) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
a matrix representation of a multiple sequence alignment, each row is a sequence, each column is a position
}
  \item{aaPrior}{
a list of two elements, "alpha"    "mix.coef", representing mixture Dirichlet prior
}
  \item{tau}{
numeric
}
  \item{priorFileName}{
string
}

}
