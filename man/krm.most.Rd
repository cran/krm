\name{krm.most}
\alias{krm.most}
\title{
Kernel-based Regression Model Maximum of adjusted Score Test
}
\description{
Performs maximum of adjusted score test for kernel-based regression models. Both Euclidean and protein sequence covariates can be used to form kernels.

}
\usage{
krm.most (formula, data, regression.type=c("logistic","linear"), 
    kern.type=c("rbf","mi","mm","prop"), n.rho=10, range.rho=0.99, 
    n.mc=2000, 
    seq.file.name=NULL, formula.kern=NULL, seq.start=NULL, seq.end=NULL,
    inference.method=c("parametric.bootstrap", "perturbation", "Davies"),
    verbose=FALSE) 

}
\arguments{
  \item{formula}{
a formula object describing the null model
}
  \item{data}{
data frame
}
  \item{regression.type}{
logistic regression or linear regression
}
  \item{kern.type}{
rbf: radial basis function kernel, a kernel type for Euclidean covariates. 
The other three kernels are for protein sequence covariates (Fong et al. 2014). mm: match-mismatch, prop: physicochemical properties, mi: mutual information, 
}
  \item{n.rho}{
integer. Number of rhos to maximize over
}
  \item{range.rho}{
numeric. A number between 0 and 1. It controls the range of rhos to use to compute kernel
}
  \item{seq.file.name}{
There are two ways to provide protein sequence information. One is to supply a sequence file named 'seq.file.name', which contains sequences in fasta format.
Two is to supply a formula through formula.kern, and the variable name should be part of data.
}
  \item{formula.kern}{
The formula for the covariates used to form the kernel. It may specify Euclidean covariates or a string covariate that contains protein sequences.
}
  \item{seq.start}{
integer. Start position of subsequence to be used in computing kernel. Only supported when the sequence is specified through formula.kern.
}
  \item{seq.end}{
integer. End position of subsequence to be used in computing kernel. Only supported when the sequence is specified through formula.kern.
}
  \item{n.mc}{
integer. Number of bootstrap samples used to compute p-values.
}
  \item{inference.method}{
parametric.bootstrap implements methods from Fong et al. (2014). 
perturbation uses methods from Wu et a. (2013)
Davies uses the upper bound method from Davies (1987) and Liu et al. (2008). 
}
  \item{verbose}{
boolean
}
}
\references{
Fong, Y. and Datta, S. and Georgiev, I. and Kwong, P. and Tomaras, G. (2014) Kernel-based
Logistic Regression Model for Protein Sequence Without Vectorialization, Biostatistics.

Wu et al. (2013) Kernel Machine SNP-Set Testing Under Multiple Candidate Kernels, Genetic epidemiology.

Davies, R. (1987) Hypothesis testing when a nuisance parameter is present only under the
alternative, Biometrika, 74, 33-43.

Liu, D., Ghosh, D. and Lin, X. (2008) Estimation and testing for the effect of a genetic
pathway on a disease outcome using logistic kernel machine regression via logistic mixed
models, BMC Bioinformatics, 9, 292.
}
\examples{

# in addition to the examples listed here, there are more examples 
# under folder R/library/krm/unitTests

\dontrun{
# the examples are not run during package build because it takes a little too long to run

# an Euclidean kernel example from Liu et al. (2008)
data=sim.liu.2008 (n=100, a=.1, seed=1) 
test = krm.most(y~x, data, formula.kern=~z.1+z.2+z.3+z.4+z.5, kern.type="rbf")


# a protein sequence kernel example
dat.file.name=paste(system.file(package="krm")[1],'/misc/y1.txt', sep="") 
seq.file.name=paste(system.file(package="krm")[1],'/misc/sim1.fasta', sep="")
dat=read.table(dat.file.name); names(dat)="y"
test = krm.most (y~1, dat, seq.file.name=seq.file.name, kern.type="mi")


} 


}
