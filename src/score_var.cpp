#include <sstream>
#include <cmath>

extern "C" {  

void score_var (double * _K, int* _n, double * _mu, double * _mom) {

    int n = *_n;
    //double (*K)[n]   = (double (*)[n])_K;	//not allowed by CRAN
    double * mu = _mu;	
    
    double mom=0;

    for (int i=0; i<n; i++) {
        mom = mom + (mu[i]*pow(1-mu[i],4) + (1-mu[i])*pow(mu[i],4)) * _K[i*n+i] * _K[i*n+i];    
    }

    for (int i=0; i<n; i++) {
        for (int k=0; k<n; k++) {
            if (k==i) continue;
            mom = mom + mu[i]*(1-mu[i]) * mu[k]*(1-mu[k]) * _K[i*n+i] * _K[k*n+k];
        }
    }

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (j==i) continue;
            mom = mom + mu[i]*(1-mu[i]) * mu[j]*(1-mu[j]) * pow(_K[i*n+j],2) * 2;
        }
    }
   
    *_mom = mom;
    
} // end 

} // end extern




