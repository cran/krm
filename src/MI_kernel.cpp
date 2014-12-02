#include <sstream>

#include "U.h"
#include "R.h"

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

void MI_kernel_str (char** seq, int* nSeq, int* seqLength, double * _tau, double * K) {

	stringstream priorStream;
	priorStream << ProteinSequenceCloud9PriorString;
    
	vector<string> seqString;
	for (int i=0; i<*nSeq; i++) {
		seqString.push_back(string(*seq+i*(*seqLength), *seqLength)); 
    }
	
    std::ostringstream oss;
    ProteinSequence * ps = new ProteinSequence (seqString, priorStream, oss);
    
    double tau = *_tau;
    
    ps->getMIKernel(K, tau);
   
    
} // end 

} // end extern




