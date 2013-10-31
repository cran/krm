\name{calcPairwiseIdentity}
\alias{calcPairwiseIdentity}
\alias{alignment2count}
\alias{alignment2trancount}
\alias{removeGap}


\title{
Functions Related to Sequence Alignment
}
\description{
Functions related to sequence alignment
}
\usage{
calcPairwiseIdentity(alignment, dissimilarity, removeGap)
alignment2count (alignment, level=20, weight=rep(1,nrow(alignment)))
alignment2trancount (alignment, weight=rep(1,nrow(alignment)))
removeGap (seq) 


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{matrix of arabic representation of sequences (1 based)}
  \item{dissimilarity}{Boolean.}
  \item{removeGap}{Boolean}
  \item{level}{integer. Size of alphabet}
  \item{weight}{numeric vector. Weights given to each sequence}

  \item{seq}{string. A string of amino acids}

}
\value{
alignment2count return T by 20 matrix, where T is the number of column in the alignment.
alignment2trancount return a T by 4 matrix, each row is the count of MM, MD, DM, DD for each position.



}
