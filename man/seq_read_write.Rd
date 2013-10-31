\name{readFastaFile}

\alias{readFastaFile}
\alias{writeFastaFile}
\alias{aa2arabic}
\alias{string2arabic}
\alias{fastaFile2arabicFile}
\alias{selexFile2arabicFile}
\alias{stringList2arabicFile}
\alias{arabic2arabicFile}
\alias{readSelexFile}
\alias{readSelexAsMatrix}
\alias{arabic2fastaFile}
\alias{readArabicFile}
\alias{readBlockFile}

\title{
Read a Fasta Sequence File
}
\description{
Read a Fasta Sequence File
}
\usage{
readFastaFile(fileName, sep = " ")
writeFastaFile (seqList, fileName)
aa2arabic (seq1) 
string2arabic (seqList) 
fastaFile2arabicFile (fastaFile, arabicFile, removeGapMajor=FALSE)
selexFile2arabicFile (selexFile, arabicFile, removeGapMajor=FALSE)
stringList2arabicFile (seqList, arabicFile, removeGapMajor=FALSE)
arabic2arabicFile (alignment, arabicFile) 
readSelexFile (fileName)
readSelexAsMatrix (fileName) 
arabic2fastaFile (alignment, fileName)
readArabicFile (fileName) 
readBlockFile (fileName) 


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{string}
  \item{fastaFile}{string}
  \item{arabicFile}{string}
  \item{selexFile}{string}
  \item{sep}{string}

  \item{seq1}{string. A string of amino acids}
  \item{seqList}{list of string.}
  \item{removeGapMajor}{Boolean}
  \item{alignment}{matrix of arabic representation of sequences (1 based)}

}
\value{
string2arabic returns a matrix of arabic numbers representing aa.
readSelexFile return a list of strings.
readArabicFile return a matrix of n by p alignment.

}
\examples{

library(RUnit)
fileName=paste(system.file(package="krm")[1],'/misc/SETpfamseed_aligned_for_testing.fasta', sep="")
seqs = readFastaFile (fileName, sep=" ")
checkEquals(length(seqs),11)

}
