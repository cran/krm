#ifndef Observable_H
#define Observable_H

#include "U.h"
#include "RandomPlus.h"

// fitnessZType is probably not the best name, 
enum FitnessZType {
  INTEGRATED, // including both finite and infinite mixture model induced prior
  CLASSIFICATION, // classificaiton likelihood, no likelihood contribution from z
  MIXTURE // mixture likelihood, observation-centric
};


class Observable {

public:
	int n; // n needs to be visible to PartitionSampler	
	int T; // dimension of each observation
	
	// return fit, e.g. integrated likelihood, for jth cluster in partition z
	virtual double getClusterFit (const int* z, int j) const = 0; 
	
	virtual double getClassificationLik (const int* z, int j) const = 0; 
	
	virtual double getMixtureLik (const int* z, int k) const = 0;
	
	virtual void writeComponentParam (const int* z, int k, ofstream & ofile) const=0;
	
	void getDistanceMatrix (double ** distanceMatrix) const {
	    for (int i=0; i<n; i++) 
			for (int j=0; j<n; j++) 
				distanceMatrix[i][j]=pairwiseDistance(i, j);
	}
	
	// fill distanceMatrix by calculating pairwise distance using a random subset of the observable dimensions
	void getRandomDistanceMatrix (double ** distanceMatrix) const {
		int numOfPositions=max(1, T/4);
		std::vector<int> positions(numOfPositions);
	    RandomPlus::SampleWithoutReplacement(T, positions, true);
	    for (int i=0; i<n; i++) 
			for (int j=0; j<n; j++) 
				distanceMatrix[i][j]=pairwiseDistance(i, j, numOfPositions, positions);
		//U::print2Dmatrix(distanceMatrix, n, n, logFile);
	}
	
protected:
	// to discourage too much jumping around or decrease the temperature in the sense of simulated annealing
	double scale;
	
	// compute pairwise distance
	virtual double pairwiseDistance (int i1, int i2) const = 0;
	virtual double pairwiseDistance (int i1, int i2, int dimension, std::vector<int> & positions) const = 0;	//changed 2013-11-01

};

#endif
