#ifndef MixtureDirichletRV_h
#define MixtureDirichletRV_h

#include "U.h"

class MixtureDirichletRV {

public: 

    MixtureDirichletRV (istream& priorStream);
    MixtureDirichletRV (string aaPriorFileName);
	MixtureDirichletRV (const MixtureDirichletRV& mrv );
	MixtureDirichletRV& operator = (const MixtureDirichletRV& mrv);
    ~MixtureDirichletRV ();

    double setWithCountsPrior(int *counts, const MixtureDirichletRV& prior); 

    void draw(double * mu, bool Log) const;

    double logDensity(double *input, bool logInput=true) const;
    double logIntegratedLik (int * y) const;
    
    double getSignature () const;
    void show(ofstream & out) const;
    void sanityCheck() const;
    
    void scaleAlpha (double tau);

private:

    int dimension;
    int mixtureOrder;
    double** alpha; // alpha[1,...,mixtureOrder][1,...,dimension]
    double *mixtureCoef;

	void readFromPriorStream (istream& source);
    void allocateMemory(int dimension, int mixtureOrder);
	void dispose();//used by destructor and by assignment operator
	void copy(const MixtureDirichletRV & rv);

};

#endif
