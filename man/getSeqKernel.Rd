\name{getSeqKernel}
\alias{getSeqKernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Protein Sequence Kernels
}
\description{
Get mutual information and other kernels for protein sequences
}
\usage{
getSeqKernel (sequences, kern.type=c("mm", "prop", "mi"), tau, call.C=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequences}{
String or list. If string, the name of a fasta file containing aligned sequences. If list, a list of strings, each string is a protein sequence. If list, call.C will be set to FALSE internally because C/C++ function needs sequence file name as input
}
  \item{kern.type}{
string. Type of kernel. mm: match-mismatch, prop: physicochemical properties, mi: mutual information.
}
  \item{tau}{
Numeric. It is the same as rho^-2.
}
  \item{call.C}{
Boolean. If TRUE, do a .C call. If FALSE, the implementation is in R. The .C call is 50 times faster.
}
}
\details{
call.C option is to allow comparison of R and C implementation. The two should give the same results and C implementation is 50 times faster. 

when kern.type is mi and call.C is TRUE and when running on linux, this function will print messages like "read ...". This message is generated from U::openRead
}
\examples{

fileName=paste(system.file(package="krm")[1],'/misc/SETpfamseed_aligned_for_testing.fasta',
     sep="")
K=getSeqKernel (fileName, kern.type="mi", tau=1, call.C=TRUE)
K


}
