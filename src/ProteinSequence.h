#ifndef ProteinSequence_H
#define ProteinSequence_H

#include "Observable.h"

#include "U.h"
#include "DirichletRV.h"
#include "MixtureDirichletRV.h"

#define C20 20

typedef int (*EMI_COUNT)[20] ;
typedef double (*EMI_PROB)[20] ;
typedef vector<MixtureDirichletRV> EMI_DIST;

typedef int (*TRAN_COUNT)[2] ;
typedef double (*TRAN_PROB)[2];
typedef vector<DirichletRV> TRAN_DIST;

// transition prior
const double taoHyperParam [2] = {.7939, .0135};
const double nuHyperParam [2] = {.9002, .5630};

// emission prior string
const string ProteinSequenceBgPriorString="  \
1 20 \
1 \
0.076 0.017 0.053 0.063 0.041 0.068 0.022 0.057 0.06 0.093 0.024 0.045 0.049 0.04 0.052 0.072 0.057 0.065 0.013 0.032 \
";

const string ProteinSequenceFlatPriorString="  \
1 20 \
1 \
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
";

const string ProteinSequenceAlpha2PriorString="  \
1 20 \
1 \
2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 \
";

const string ProteinSequenceAlpha20PriorString="  \
1 20 \
1 \
20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 \
";

const string ProteinSequenceCloud9PriorString="  \
9 20 \
0.178091 \
0.270671 0.039848 0.017576 0.016415 0.014268 0.131916 0.012391 0.022599 0.020358 0.030727 0.015315 0.048298 0.053803 0.020662 0.023612 0.216147 0.147226 0.065438 0.003758 0.009621 \
0.056591 \
0.021465 0.0103 0.011741 0.010883 0.385651 0.016416 0.076196 0.035329 0.013921 0.093517 0.022034 0.028593 0.013086 0.023011 0.018866 0.029156 0.018153 0.0361 0.07177 0.419641 \
0.0960191 \
0.561459 0.045448 0.438366 0.764167 0.087364 0.259114 0.21494 0.145928 0.762204 0.24732 0.118662 0.441564 0.174822 0.53084 0.465529 0.583402 0.445586 0.22705 0.02951 0.12109 \
0.0781233 \
0.070143 0.01114 0.019479 0.094657 0.013162 0.048038 0.077 0.032939 0.576639 0.072293 0.02824 0.080372 0.037661 0.185037 0.506783 0.073732 0.071587 0.042532 0.011254 0.028723 \
0.0834977 \
0.041103 0.014794 0.00561 0.010216 0.153602 0.007797 0.007175 0.299635 0.010849 0.999446 0.210189 0.006127 0.013021 0.019798 0.014509 0.012049 0.035799 0.180085 0.012744 0.026466 \
0.0904123 \
0.115607 0.037381 0.012414 0.018179 0.051778 0.017255 0.004911 0.796882 0.017074 0.285858 0.075811 0.014548 0.015092 0.011382 0.012696 0.027535 0.088333 0.94434 0.004373 0.016741 \
0.114468 \
0.093461 0.004737 0.387252 0.347841 0.010822 0.105877 0.049776 0.014963 0.094276 0.027761 0.01004 0.187869 0.050018 0.110039 0.038668 0.119471 0.065802 0.02543 0.003215 0.018742 \
0.0682132 \
0.452171 0.114613 0.06246 0.115702 0.284246 0.140204 0.100358 0.55023 0.143995 0.700649 0.27658 0.118569 0.09747 0.126673 0.143634 0.278983 0.358482 0.66175 0.061533 0.199373 \
0.234585 \
0.005193 0.004039 0.006722 0.006121 0.003468 0.016931 0.003647 0.002184 0.005019 0.00599 0.001473 0.004158 0.009055 0.00363 0.006583 0.003172 0.00369 0.002967 0.002772 0.002686 \
";

class ProteinSequence : public Observable {

public: 
	ProteinSequence (string sequenceFileName, stringstream& priorStream, ostream& _logFile );
	ProteinSequence (vector<string> seqString, stringstream& priorStream, ostream& _logFile );
	~ProteinSequence ();
	
	double getClusterFit (const int* z, int j) const;
	double getMixtureLik (const int* z, int k) const;
	double getClassificationLik (const int* z, int k) const;
	void writeComponentParam (const int* z, int k, ofstream & ofile) const;
	void getMIKernel (double * K, double tau); // muPrior is changed
    double hmmMargLlik (int i, int j) const;
	
private:
	int ** sequence;
	ostream& logFile;

    MixtureDirichletRV muPrior;
    DirichletRV taoPrior;
    DirichletRV nuPrior;

	void readSeqNamesFromFasta (const char * fileName, vector<string>& seqNames) const;
	void readFastaFile (string fileName);
	void setSequence (vector<string> seqString);
	int getAAcount(const int *z, short j, int * count) const;
    void getAAcount(int i, int j, int * count) const;
	void getTranCount(const int *z, short j, int * taoCount, int * nuCount) const;
    void getTranCount(int i, int j, int * taoCount, int * nuCount) const;

	double pairwiseDistance (int i1, int i2) const;
	double pairwiseDistance (int i1, int i2, int length, std::vector<int> & positions) const;
};

#endif
