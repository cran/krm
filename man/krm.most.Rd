\name{krm.most}
\alias{krm.most}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel-based Regression Model Maximum of adjusted Score Test
}
\description{
Computes maximum of adjusted score test. Obtain p value through parametric bootstrap

}
\usage{
krm.most (formula, data, regression.type=c("logistic","linear"), 
    kern.type=c("rbf","mi","mm","prop"), 
    n.rho=10, range.rho=0.99, n.mc=2000, seq.file.name=NULL, formula.kern=NULL,
    inference.method=c("parametric.bootstrap", "perturbation", "LGL2008"),
    verbose=FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
a formula object describing the null model
}
  \item{data}{
data frame
}
  \item{regression.type}{
a string
}
  \item{formula.kern}{
formula. The formula for the covariates used to form the kernel
}
  \item{seq.file.name}{
string. Name of a file containing sequences in fasta format
}
  \item{kern.type}{
string. Type of kernel. mm: match-mismatch, prop: physicochemical properties, mi: mutual information, rbf: radial basis function
}
  \item{n.rho}{
integer. Number of rhos to maximize over
}
  \item{range.rho}{
numeric. A number between 0 and 1. It controls the range of rhos to use to compute kernel
}
  \item{n.mc}{
integer. Number of bootstrap samples
}
  \item{inference.method}{
string
}
  \item{verbose}{
boolean
}
}
\examples{

\dontrun{
# the examples are not run during package build because it takes a little too long to run

data=sim.liu.2008 (n=100, a=.1, seed=1) 
test = krm.most(y~x, data, formula.kern=~z.1+z.2+z.3+z.4+z.5, kern.type="rbf")


dat.file.name=paste(system.file(package="krm")[1],'/misc/y1.txt', sep="") 
seq.file.name=paste(system.file(package="krm")[1],'/misc/sim1.fasta', sep="")
dat=read.table(dat.file.name); names(dat)="y"
test = krm.most (y~1, dat, seq.file.name=seq.file.name, kern.type="mi")


} 


}
