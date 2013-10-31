#include "DirichletRV.h"
#include "U.h"
#include "RandomPlus.h"
#include "MathPlus.h"
#include "StatPlus.h"

DirichletRV::DirichletRV (short _dimension){
	set(_dimension, defaultAlpha);
}

DirichletRV::DirichletRV (short _dimension, double _alpha){
	set(_dimension, _alpha);
}

DirichletRV::DirichletRV (short _dimension, const double * _alpha){
	set(_dimension, _alpha);
}

DirichletRV::DirichletRV (const DirichletRV& rv) {
	set(rv.dimension, rv.alpha);
}

DirichletRV& DirichletRV::operator = (const DirichletRV& rv) {
    delete[] alpha;
	set(rv.dimension, rv.alpha);
	return *this;	
}

DirichletRV::~DirichletRV () {
    delete[] alpha;
}

void DirichletRV::set(short _dimension, double _alpha0) {
	std::vector<double> _alpha(_dimension);
    for (int i=0; i<_dimension; i++) _alpha[i]=_alpha0;
	set(_alpha);
}

void DirichletRV::set(std::vector<double>& _alpha) {
	dimension=_alpha.size();
    alpha=new double[dimension];
    for (int i=0; i<dimension; i++) alpha[i]=_alpha[i];	
}

void DirichletRV::set(short _dimension, const double* _alpha) {
	dimension=_dimension;
    alpha=new double[dimension];
    for (int i=0; i<dimension; i++) alpha[i]=_alpha[i];	
}

// this works under the assumption that dimension is unchanged
void DirichletRV::setWithCountsPrior(int *counts, const DirichletRV & prior) {
    for (int i=0; i<dimension; i++) alpha[i]=prior.alpha[i] + counts[i];
}

void DirichletRV::draw(double * mu, bool Log){
    RandomPlus::rdirichlet(alpha, dimension, mu, Log);
}

double DirichletRV::logDensity(double * input, bool logInput) const {
    if (logInput) return StatPlus::ddirichlet_log(input, alpha, dimension, true);
    else return StatPlus::ddirichlet(input, alpha, dimension, true);
}

double DirichletRV::logIntegratedLik (int * y) const {
	std::vector<double> newAlpha(dimension);
	for (int m=0; m<dimension; m++) newAlpha[m] = alpha[m] + y[m];
	return MathPlus::lbeta(newAlpha) - MathPlus::lbeta(alpha, dimension);
}

void DirichletRV::show() const{
//	cout << "BEGIN DirichletRV::show\n";
	//U::printArray(alpha, dimension); // comment out to resolve R CMD check complaint
//	cout << "END DirichletRV::show\n";
}
