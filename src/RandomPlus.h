#ifndef RandomPlus_H
#define RandomPlus_H

#include "U.h"
#include "rvgs.h"

class RandomPlus
{

public:
	
    // from 1:n pick k numbers and write it to y 
	template <class T>
    static void SampleWithoutReplacement(int n, int k, T* y, bool zeroBased=false) // Change from T n, T k to int n, int k (2013-11-01)
    {
        T i, j;
        T x[n]; // a temp array
        for (i = 0; i < n; i++)
    	    x[i] = i;
        for (i = 0; i < k; i++) {
        	j = (T) (n * myunif_rand());
        	if (!zeroBased) 
				y[i] = x[j] + 1;
        	else
				y[i] = x[j];
        	x[j] = x[--n];
        }
    }
	template <class T>
    static void SampleWithoutReplacement(int n, std::vector<T> & y, bool zeroBased=false)
    {
        int k = y.size();
        T i, j;
        std::vector<T> x(n); // a temp array
        for (i = 0; i < n; i++)
    	    x[i] = i;
        for (i = 0; i < k; i++) {
        	j = (T) (n * myunif_rand());
        	if (!zeroBased) 
				y[i] = x[j] + 1;
        	else
				y[i] = x[j];
        	x[j] = x[--n];
        }
    }

	static void rdirichlet(double* alpha, short categories, double* out, bool Log){
		for (int i=0; i<categories; i++) out[i] = rGamma(alpha[i], Log);
		if (Log) {
			double logSum=U::logSumExp(out, categories);
			for (int i=0; i<categories; i++) out[i] -= logSum;
		} else {
			double sum=U::sum(out, categories);
			for (int i=0; i<categories; i++) out[i] /= sum;
		}
	}

	// return 0, 1, ... , dimension-1	
	static int rGeneralizedBern(double * mean, int dimension) {

		if (dimension==1) return 0;
		double p=myunif_rand();
		double cusum=0;
		int d=0;
		for (; d<dimension; d++){
			cusum+=mean[d];
			if (cusum>p) break;
		}
		return d;
	}
	static int rGeneralizedBern(vector<double>& mean) {

		int dimension=mean.size();

		if (dimension==1) return 0;
		double p=myunif_rand();
		double cusum=0;
		int d=0;
		for (; d<dimension; d++){
			cusum+=mean[d];
			if (cusum>p) break;
		}
		return d;
	}
    // normalize logMean and draw from a generalized Bernoulli
	static int drawGeneralizedBern (double * logMean0, int numCells){
	    double logSumOfMean = U::logSumExp(logMean0, numCells);
        std::vector<double> mean(numCells);
	    for (int j=0; j<numCells; j++) mean[j] = exp(logMean0[j]-logSumOfMean);
	    return RandomPlus::rGeneralizedBern(mean);
	}
	// the only difference from last function is that it assigns a value to logProb
	static int drawGeneralizedBern (double * logMean0, int numCells, double & logProb){
	    double logSumOfMean = U::logSumExp(logMean0, numCells);
	    std::vector<double> mean(numCells);
	    for (int j=0; j<numCells; j++) mean[j] = exp(logMean0[j]-logSumOfMean);
	    int draw = RandomPlus::rGeneralizedBern(mean);
	    logProb = logMean0[draw]-logSumOfMean;
	    return draw;
	}

};    
#endif
