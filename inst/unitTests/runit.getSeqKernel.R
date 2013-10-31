### --- Test setup ---

if(FALSE) {
  library("RUnit")
  library("krm")
}

test.getSeqKernel <- function() {


tolerance=1e-3
# more stringent tolerance for one system to ensure algorithm accuracy
if (R.Version()$system %in% c("x86_64, mingw32")) {
    tolerance=1e-6
} 
RNGkind("Mersenne-Twister", "Inversion")


fileName=paste(system.file(package="krm")[1],'/misc/SETpfamseed_aligned_for_testing.fasta', sep="")

K=getSeqKernel (fileName, kern.type="mi", tau=.01, call.C=T)
checkEqualsNumeric(
    c(K[1:2, 1:2])
    , 
    c(1.0000000, 0.1038921, 0.1038921, 1.0000000)
    , tolerance = tolerance
)

## compare R and C implementation, commented out b/c it takes a while to run
#K.1=getSeqKernel (fileName, kern.type="mi", tau=.01, call.C=FALSE)
#checkEqualsNumeric(K, K.1, tolerance = tolerance)

K=getSeqKernel (fileName, kern.type="mi", tau=1, call.C=T)
checkEqualsNumeric(
    c((K^0.01)[1:2, 1:2])
    , 
    c(1.0000000, 0.1038921, 0.1038921, 1.0000000)
    , tolerance = tolerance
)

seq.alignment <- readFastaFile(fileName)
K=getSeqKernel (seq.alignment[1:2], kern.type="mi", tau=.01)
checkEqualsNumeric(
    c(K[1:2, 1:2])
    , 
    c(1.0000000, 0.1038921, 0.1038921, 1.0000000)
    , tolerance = tolerance
)


aaList=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
sequences=as.list(aaList[-21])
names(sequences)=aaList[-21]
K.aa=getSeqKernel (sequences, kern.type="mi", tau=1)


}
