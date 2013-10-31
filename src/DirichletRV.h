#ifndef DirichletRV_h
#define DirichletRV_h

#include "U.h"

const int defaultAlpha=1;

class DirichletRV {

public: 

    DirichletRV (short _dimension);
    DirichletRV (short _dimension, double _alpha);
    DirichletRV (short _dimension, const double * _alpha);
    DirichletRV (const DirichletRV& rv);
	DirichletRV& operator = (const DirichletRV& rv);
    ~DirichletRV ();
    
	// change the distribution
    void setWithCountsPrior(int *counts, const DirichletRV & prior); 

    void draw(double * mu, bool Log);
    
    double logDensity(double * input, bool logInput) const;  
    double logIntegratedLik (int * y) const;
    
    void show () const;

private:
	
    short dimension;
    double* alpha; 
    
    void set(short _dimension, const double * _alpha);
    void set(short _dimension, double _alpha);
    void set(std::vector<double>& _alpha);
};

#endif
