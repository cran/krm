#include "MixtureDirichletRV.h"
#include "U.h"
#include "RandomPlus.h"
#include "StatPlus.h"
#include "MathPlus.h"

MixtureDirichletRV::MixtureDirichletRV (istream& priorStream) {
	readFromPriorStream(priorStream);
}

MixtureDirichletRV::MixtureDirichletRV (string aaPriorFileName) {
    ifstream aaPriorFile;
    std::ostringstream oss;
    U::openRead (aaPriorFile, aaPriorFileName, oss);
	readFromPriorStream(aaPriorFile);
    aaPriorFile.close();
}   

MixtureDirichletRV::MixtureDirichletRV (const MixtureDirichletRV& rv ) {
    allocateMemory (rv.dimension, rv.mixtureOrder);
    copy(rv);
}

MixtureDirichletRV& MixtureDirichletRV::operator = (const MixtureDirichletRV& rv) {
    if(rv.dimension!=dimension || rv.mixtureOrder!=mixtureOrder) dispose();
	copy(rv);
	return *this;	
}

MixtureDirichletRV::~MixtureDirichletRV () {
	dispose();
}

void MixtureDirichletRV::dispose(){
    delete [] mixtureCoef;
    for (int d=0; d<mixtureOrder; d++) delete [] alpha[d];
    delete [] alpha;    
}

void MixtureDirichletRV::readFromPriorStream (istream& source) {
	int _mixtureOrder, _dimension;
    source >> _mixtureOrder >> _dimension;
    allocateMemory (_dimension, _mixtureOrder);
            
    for (int d=0; d<mixtureOrder; d++) {
        source >> mixtureCoef[d];
        for (int m=0; m<dimension; m++) source >> alpha[d][m];     
    }       
    #ifdef __DEBUG__
    //U::print2Dmatrix(alpha, mixtureOrder, dimension);
    #endif  
}

void MixtureDirichletRV::allocateMemory(int _dimension, int _mixtureOrder) {
	dimension=_dimension;
	mixtureOrder=_mixtureOrder;
    mixtureCoef = new double[mixtureOrder];
    alpha = new double*[mixtureOrder];
    for (int d=0; d<mixtureOrder; d++) alpha[d] = new double[dimension];
}

void MixtureDirichletRV::copy(const MixtureDirichletRV & rv) {
    for (int d=0; d<mixtureOrder; d++) {
        mixtureCoef[d]=rv.mixtureCoef[d];
        for (int m=0; m<dimension; m++) alpha[d][m]=rv.alpha[d][m];     
    }       
}

/*
void MixtureDirichletRV::expandMixture (double proportionOfExisting, double _alpha){
	int newMixtureOrder=mixtureOrder + 1;
	double * newMixtureCoef = new double[newMixtureOrder];
	for (int i=0; i<mixtureOrder; i++) newMixtureCoef[i] = mixtureCoef[i]*proportionOfExisting;
	for (int i=0; i<1; i++) 
		newMixtureCoef[mixtureOrder+i] = 1-proportionOfExisting;
	
	double ** newAlpha=new double*[newMixtureOrder];
	for (int i=0; i< mixtureOrder; i++) newAlpha[i] = alpha[i];
	for (int i=0; i< 1; i++) {
		newAlpha[mixtureOrder+i] = new double[dimension];
		for (int m=0; m<dimension; m++) newAlpha[mixtureOrder+i][m]=_alpha;
	}

	mixtureOrder = newMixtureOrder;
	delete[] mixtureCoef;
	mixtureCoef=newMixtureCoef;
	delete[] alpha;
	alpha=newAlpha;
	
	//U::printArray(mixtureCoef, mixtureOrder);
	//U::print2Dmatrix(alpha, mixtureOrder, dimension);	
}

void MixtureDirichletRV::expandMixture (double proportionOfExisting, MixtureDirichletRV& toAdd){
	int newMixtureOrder=mixtureOrder + toAdd.mixtureOrder;
	double * newMixtureCoef = new double[newMixtureOrder];
	for (int i=0; i<mixtureOrder; i++) newMixtureCoef[i] = mixtureCoef[i]*proportionOfExisting;
	for (int i=0; i<toAdd.mixtureOrder; i++) 
		newMixtureCoef[mixtureOrder+i] = toAdd.mixtureCoef[i]*(1-proportionOfExisting);
	
	double ** newAlpha=new double*[newMixtureOrder];
	for (int i=0; i< mixtureOrder; i++) newAlpha[i] = alpha[i];
	for (int i=0; i< toAdd.mixtureOrder; i++) {
		newAlpha[mixtureOrder+i] = new double[dimension];
		U::copyArray(toAdd.alpha[i], newAlpha[mixtureOrder+i], dimension);
	}

	mixtureOrder = newMixtureOrder;
	delete[] mixtureCoef;
	mixtureCoef=newMixtureCoef;
	delete[] alpha;
	alpha=newAlpha;
	
	//U::printArray(mixtureCoef, mixtureOrder);
	//U::print2Dmatrix(alpha, mixtureOrder, dimension);
}
*/

void MixtureDirichletRV::draw(double * mu, bool Log) const {
	int d=RandomPlus::rGeneralizedBern(mixtureCoef, mixtureOrder);
    RandomPlus::rdirichlet(alpha[d], dimension, mu, Log);
}

// this function works under the assumption that there is no need to allocate memory
// returns integrated likelihood as a byproduct
double MixtureDirichletRV::setWithCountsPrior(int *counts, const MixtureDirichletRV& prior) {
    for (int d=0; d<mixtureOrder; d++) {
        for (int m=0; m<dimension; m++) alpha[d][m]=prior.alpha[d][m]+counts[m];           
        mixtureCoef[d]=exp( log(prior.mixtureCoef[d]) +
                            MathPlus::lbeta(alpha[d], (int)dimension) -
                            MathPlus::lbeta(prior.alpha[d], (int)dimension) 
                        );
    }
    double marginalLik = U::sum(mixtureCoef, mixtureOrder);
    for (int d=0; d<mixtureOrder; d++) mixtureCoef[d]/=marginalLik;
    return marginalLik;
}

double MixtureDirichletRV::logDensity(double * input, bool logInput) const {
    std::vector<double> logDensities(mixtureOrder);  
    if (logInput)
	    for (int d=0; d<mixtureOrder; d++) 
			logDensities[d]=StatPlus::ddirichlet_log(input, alpha[d], dimension, true);
	else
	    for (int d=0; d<mixtureOrder; d++) 
			logDensities[d]=StatPlus::ddirichlet(input, alpha[d], dimension, true);
    return U::logSumExp(logDensities, mixtureCoef);
}

double MixtureDirichletRV::logIntegratedLik (int * y) const {
	std::vector<double> toAdd(mixtureOrder);
	std::vector<double> newAlpha(dimension);
    for (int d=0; d<mixtureOrder; d++) {
	    for (int m=0; m<dimension; m++) newAlpha[m] = alpha[d][m] + y[m];
		toAdd[d] = MathPlus::lbeta(newAlpha) - MathPlus::lbeta(alpha[d], dimension);
	}
	return U::logSumExp (toAdd, mixtureCoef);
}

// return the sum of all elements of alpha, for diagnostic purpose
double MixtureDirichletRV::getSignature() const {
    double out=0;
    for (int d=0; d<mixtureOrder; d++) 
        for (int m=0; m<dimension; m++) 
            out += pow((-1.0),m)*alpha[d][m];
    out += mixtureCoef[0];
    return out;
}

void MixtureDirichletRV::show(ofstream& out) const {
	out << "showing MixtureDirichletRV ...\n";
    for (int d=0; d<mixtureOrder; d++) {
		out << mixtureCoef[d] << endl;
        for (int m=0; m<dimension; m++) out << alpha[d][m] << " ";
        out << endl;
	}
}

void MixtureDirichletRV::sanityCheck() const {
	//std::vector<double> a(dimension);
	//for (int i=0; i<dimension; i++) a[i]=log(1.0/dimension);
	//cout << "log density of the diagnoal vector is: " << logDensity(a, true);// commented out b/c R CMD check complains
}


void MixtureDirichletRV::scaleAlpha (double tau) {
    for (int d=0; d<mixtureOrder; d++) {
	    for (int m=0; m<dimension; m++) 
            alpha[d][m] = alpha[d][m] * tau;
	}     
}
