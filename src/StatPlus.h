#ifndef StatPlus_H
#define StatPlus_H

#include "MathPlus.h"

class StatPlus {
public:
	static double loggedNormalDensity (double x, double sigma) {
		return  - .5*log (2*_PI) - log(sigma) -.5*pow(x/sigma,2) ;
	}

	static double loggedNormalDensity (double x, double mu, double tao) {
		return  - .5*log (2*_PI) + .5*log(tao) -.5*pow(x-mu,2)*tao ;
	}

	static double loggedExponentialDensity (double x, double rate) {
		return  log(rate) - rate * x ;
	}

	static double loggedGammaDensity (double x, double shape, double rate) {
		double temp=0; //log( Gamma(shape) )
		for (int i=1; i<shape; i++) temp += log((double)i);		
		return (shape-1)*log(x) - rate*x - temp + shape*log(rate);
	}
	static double ddirichlet (double* x, double * alpha, short categories,  bool logOutput=false){
		double out=0;
		double sum=0;	
		for (int i=0; i<categories; i++) sum+=alpha[i];
		out+= MathPlus::mylgamma(sum);	 
		for (int i=0; i<categories; i++) out-=MathPlus::mylgamma(alpha[i]); 
		for (int i=0; i<categories; i++) {
			out+=(alpha[i]-1)*log(x[i]);
		}
		return logOutput?out:exp(out);
	}
	static double ddirichlet_log (double* logx, double * alpha, short categories,  bool logOutput=false){
		double out=0;
		double sum=0;
		for (int i=0; i<categories; i++) sum+=alpha[i];
		out+=MathPlus::mylgamma(sum);	 
		for (int i=0; i<categories; i++) out-=MathPlus::mylgamma(alpha[i]); 
		for (int i=0; i<categories; i++) {
			out+=(alpha[i]-1)*logx[i];
		}
		return logOutput?out:exp(out);
	}


		
};

#endif
