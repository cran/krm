#include <sstream>
#include <math.h>

extern "C" {  

void var_Q (double * _W, int* _n, double * _variance, double * _extra_kurtosis, double * _mom) {

    int n = *_n;
    //double (*W)[n]   = (double (*)[n])_W;	//not allowed by CRAN
    double * variance = _variance;	
    double * extra_kurtosis = _extra_kurtosis;	
    
    double mom=0;

    for (int i=0; i<n; i++) {
        mom = mom + pow(_W[i*n+i],2) * extra_kurtosis[i];
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            mom = mom + 2 * pow(_W[i*n+j],2) * variance[i] * variance[j];
        }
    }
   
    *_mom = mom;
    
} // end 

} // end extern




