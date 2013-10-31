#include <map>

#include "ProteinSequence.h"
#include "RandomPlus.h"

#define C2 2

ProteinSequence::ProteinSequence(string sequenceFileName, stringstream& priorStream, ostream& _logFile) 
: logFile(_logFile), muPrior(priorStream), taoPrior(C2, taoHyperParam), nuPrior(C2, nuHyperParam)
{
	readFastaFile (sequenceFileName);
	//muPrior.show(logFile);
}

ProteinSequence::~ProteinSequence() {
    for (int i=0; i<n; i++) 
        delete[] sequence[i];
    delete[] sequence;
}

void ProteinSequence::readFastaFile (string fileName) {
	ifstream msaFile; //multiple sequence alignment file
    U::openRead (msaFile, fileName, logFile);
	vector<string> seqString;
	string seqName, seqPart, seq;
	getline(msaFile, seqName);
	while (!msaFile.eof()) {
		seq = "";
		getline(msaFile, seqPart);
		while (!msaFile.eof() && seqPart.compare (0, 1, ">")!=0) {
			seq += seqPart;
			getline(msaFile, seqPart);
		}
		seqString.push_back(seq); 
	}
	msaFile.close();
	
	n=seqString.size();
	T=seqString[0].length();

	// map from ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-","*") to 0-20
	map<char,int> aaMap;
	int counter=0;
	aaMap['A']=counter++; aaMap['C']=counter++; aaMap['D']=counter++; aaMap['E']=counter++;
	aaMap['F']=counter++; aaMap['G']=counter++; aaMap['H']=counter++; aaMap['I']=counter++;
	aaMap['K']=counter++; aaMap['L']=counter++; aaMap['M']=counter++; aaMap['N']=counter++;
	aaMap['P']=counter++; aaMap['Q']=counter++; aaMap['R']=counter++; aaMap['S']=counter++;
	aaMap['T']=counter++; aaMap['V']=counter++; aaMap['W']=counter++; aaMap['Y']=counter++;
	// all gaps are mapped to 20
	aaMap['-']=counter; aaMap['*']=counter; aaMap['.']=counter; 

    sequence = new int*[n];
    for (int i=0; i<n; i++) 
        sequence[i] = new int[T];

	for (int i=0; i<n; i++) 
        for (int t=0; t<T; t++)
            sequence[i][t] = aaMap[toupper(seqString[i][t])];
    #ifdef __DEBUG__
	//U::print2Dmatrix(sequence, n, T, logFile);
	#endif
}

void ProteinSequence::readSeqNamesFromFasta (const char * fileName, vector<string>& seqNames) const {
	string seqName, seq;
	ifstream msaFile; //multiple sequence alignment file

	msaFile.open(fileName, ios::in);

	getline(msaFile, seqName);
	while (!msaFile.eof()) {
		getline(msaFile, seq);
		while (!msaFile.eof() && seq.compare (0, 1, ">")!=0) {
			getline(msaFile, seq);
		}
		seqNames.push_back(seqName.substr(1,seqName.find_first_of(" ")-1));
		//cout << seq; // for debug
		seqName=seq;
	}

	msaFile.close();
}

// for a given z vector, count aa for the  j'th cluster and store in count
// return size of cluster j
int ProteinSequence::getAAcount(const int *z, short j, int * count) const {
	int n_j=0;
	for (int t=0; t<T; t++) 
		for (int m=0; m<C20; m++) 
			count [t*C20+m] =0;

	for (int i=0; i<n; i++) 
		if (z[i]==j) {
			n_j ++;
			for (int t=0; t<T; t++) 
				if (sequence[i][t]!=C20)
					count[t*C20+sequence[i][t] ] ++;
		}
		
	return n_j;
}

// for two sequenced indexed by i and j, count aa 
void ProteinSequence::getAAcount(int i, int j, int * count) const {
	for (int t=0; t<T; t++) 
		for (int m=0; m<C20; m++) 
			count [t*C20+m] =0;

	for (int t=0; t<T; t++) 
		if (sequence[i][t]!=C20)
			count[t*C20+sequence[i][t] ] ++;

	for (int t=0; t<T; t++) {
		if (sequence[j][t]!=C20)
			count[t*C20+sequence[j][t] ] ++;
   }
}

// for two sequences indexed by i and j, count transition from Match 
void ProteinSequence::getTranCount(int i, int j, int * taoCount, int * nuCount) const {
	for (int t=0; t<T; t++) 
		for (int m=0; m<C2; m++) {
			taoCount [t*C2+m] =0;
			nuCount [t*C2+m] =0;
		}
		
	// use v loop to do i,j	
	for (int v=0; v<2; v++) {
        
        if (v==1) i=j; // do j
         
		// first position gets special treatment
		int t=0;
		if (sequence[i][t]!=20) 
			taoCount[t*C2+0] ++;
		else 
			taoCount[t*C2+1] ++;		
		// continue with the rest	
		for (t=1; t<T; t++) {
			if (sequence[i][t-1]!=20) {
				if (sequence[i][t]!=20) taoCount[t*C2+0] ++;
				else taoCount[t*C2+1] ++;			
			} else {
				if (sequence[i][t]!=20) nuCount[t*C2+0] ++;
				else nuCount[t*C2+1] ++;			
			}
		}
		
	}
}

// for a given z vector, count transition from Match for the j'th cluster 
void ProteinSequence::getTranCount(const int *z, short j, int * taoCount, int * nuCount) const {
	for (int t=0; t<T; t++) 
		for (int m=0; m<C2; m++) {
			taoCount [t*C2+m] =0;
			nuCount [t*C2+m] =0;
		}
	for (int i=0; i<n; i++) 
		if (z[i]==j) {
			// first position gets special treatment
			int t=0;
			if (sequence[i][t]!=20) 
				taoCount[t*C2+0] ++;
			else 
				taoCount[t*C2+1] ++;		
			// continue with the rest	
			for (t=1; t<T; t++) {
				if (sequence[i][t-1]!=20) {
					if (sequence[i][t]!=20) taoCount[t*C2+0] ++;
					else taoCount[t*C2+1] ++;			
				} else {
					if (sequence[i][t]!=20) nuCount[t*C2+0] ++;
					else nuCount[t*C2+1] ++;			
				}
			}
		}
}

double ProteinSequence::getMixtureLik (const int* z, int k) const {
	logFile << "kmeans2::getMixtureLik not implmented" << endl;
	return -1;
	//exit(-1); // commented out b/c R CMD check complains
}

double ProteinSequence::getClassificationLik (const int* z, int k) const {
	logFile << "ProteinSequence::getClassificationLik not implmented" << endl;
	return -1;
	//exit(-1); // commented out b/c R CMD check complains
}

void ProteinSequence::writeComponentParam (const int* z, int k, ofstream & ofile) const {
}

double ProteinSequence::getClusterFit (const int* z, int j) const {
	double ll=0;	
//	int aaCount[T][C20];	
//	int taoCount[T][C2];
//	int nuCount[T][C2];
   //CRAN does not allow variable length array
   int * aaCount = new int [T*C20];
   int * taoCount = new int [T*C2];
   int * nuCount = new int [T*C2];
	int n_tao, n_nu, n_aa;

	int n_j = getAAcount(z, j, aaCount);
	if (n_j==0) return 0;
	getTranCount(z, j, taoCount, nuCount); 
    for (int t=0; t<T; t++) {
		// we could allow each cluster, column to use its own prior, but when j=1000, this fails
		// thus we take the easy way out and use taoPrior[0][0] for everyone
		
		// transition likelihoods
		n_tao=0;
        for (int m=0; m<C2; m++) 
			n_tao += taoCount[t*C2+m];
        if (n_tao>0) 
			ll += taoPrior.logIntegratedLik(taoCount+t*C2);

		n_nu=0;
        for (int m=0; m<C2; m++) 
			n_nu += nuCount[t*C2+m];
        if (n_nu>0) 
			ll += nuPrior.logIntegratedLik(nuCount+t*C2);

        // emission likelihood
        n_aa=0;
        for (int m=0; m<C20; m++) 
			n_aa += aaCount[t*C20+m];
        if (n_aa>0) 
			ll += muPrior.logIntegratedLik(aaCount+t*C20);
	}
    
    delete [] aaCount;
    delete [] taoCount;
    delete [] nuCount;
	return ll;
}


// get the HMM marginal log likelihood for two sequences i and j
double ProteinSequence::hmmMargLlik (int i, int j) const {

   int * aaCount = new int [T*C20];
   int * taoCount = new int [T*C2];
   int * nuCount = new int [T*C2];
	int n_tao, n_nu, n_aa;

    double ll=0;
    
        
	getAAcount(i, j, aaCount);
	getTranCount(i, j, taoCount, nuCount); 
    for (int t=0; t<T; t++) {        		
        // emission likelihood
        n_aa=0;
        for (int m=0; m<C20; m++) 
			n_aa += aaCount[t*C20+m];
        if (n_aa>0) 
			ll += muPrior.logIntegratedLik(aaCount+t*C20);

		// transition likelihoods
		n_tao=0;
        for (int m=0; m<C2; m++) 
			n_tao += taoCount[t*C2+m];
        if (n_tao>0) 
			ll += taoPrior.logIntegratedLik(taoCount+t*C2);

		n_nu=0;
        for (int m=0; m<C2; m++) 
			n_nu += nuCount[t*C2+m];
        if (n_nu>0) 
			ll += nuPrior.logIntegratedLik(nuCount+t*C2);

    }
    
    delete [] aaCount;
    delete [] taoCount;
    delete [] nuCount;
    return ll;
     
}

void ProteinSequence::getMIKernel (double * _K, double tau) {
    
    //muPrior.scaleAlpha(tau);
     
    for (int i=0; i<n; i++) {
		for (int j=i+1; j<n; j++) {
			_K[i*n+j]=exp(tau * (hmmMargLlik(i,j) - (hmmMargLlik(i,i)+hmmMargLlik(j,j))/2) );
        }     
   }
        
    for (int i=0; i<n; i++) 
		for (int j=0; j<i; j++) 
		    _K[i*n+j]=_K[j*n+i];
		
    for (int i=0; i<n; i++) 
        _K[i*n+i]=1;
		
}

double ProteinSequence::pairwiseDistance (int i1, int i2) const {
	// countGaps is true when you want to count gaps
	// whether we count a position which is both gaps in the denominator matters for the tree topology in bottomup clustering
	// it seems that removing gaps when counting denominator is better
	bool countGaps=false;
	
	int *seq1=sequence[i1];
	int *seq2=sequence[i2];

	int diff=0;
	int actual_length=0;
	for (int t=0; t<T; t++) {
		if ((seq1[t]!=20) | (seq2[t]!=20)) actual_length++;
		if (seq1[t]!=seq2[t]) diff++;
	}
	
	if (!countGaps) {
		if (actual_length==0) {
			//cout << "LatentClass::pairwiseDistance(): actual length is 0.";  // commented out b/c R CMD check complains
			return 1;
		} else
			return (double)diff/actual_length;
	} else
		return (double)diff/T;
}

double ProteinSequence::pairwiseDistance (int i1, int i2, int length, std::vector<int> & positions) const {
	bool countGaps=false;

	int *seq1=sequence[i1];
	int *seq2=sequence[i2];

	int diff=0;
	int actual_length=0;
	for (int t=0; t<length; t++) {
		if ((seq1[positions[t]]!=20) | (seq2[positions[t]]!=20)) actual_length++;
		if (seq1[positions[t]]!=seq2[positions[t]]) diff++;
	}
	
	if (!countGaps) {
		if (actual_length==0) {
			//cout << "LatentClass::pairwiseDistance(): actual length is 0.";  // commented out b/c R CMD check complains
			return 1;
		} else
			return (double)diff/actual_length;
	} else
		return (double)diff/length;
}

