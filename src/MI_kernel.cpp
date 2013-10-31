#include <sstream>

#include "U.h"

#include "ProteinSequence.h"
#include "MathPlus.h"

double myunif_rand() { return 1; } // this is because we have included other code that needs myunif_rand() 

extern "C" {  

// this function can be tested in cbclust main

void MI_kernel (char** dataFileCh, int* dataFileLen, double * _tau, double * K) {

    string dataFile=string(*dataFileCh, *dataFileLen);

	stringstream priorStream;
	priorStream << ProteinSequenceCloud9PriorString;
    
    std::ostringstream oss;
    ProteinSequence * ps = new ProteinSequence (dataFile, priorStream, oss);
    
    double tau = *_tau;
    
    ps->getMIKernel(K, tau);
   
    
} // end 

} // end extern




