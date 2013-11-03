#ifndef MathPlus_H
#define MathPlus_H

#include "math.h"

class MathPlus {
public:
	static double mylgamma(double a) {
		if (a==0) return 0;
		#ifdef __LINUX__
		int sign[1];
		return lgamma_r(a, sign);
		#else
		return lgamma(a);
		#endif		
	}	
	
	static double lbeta (double* alpha, int m, bool logged=true) {
		double out=0;
		double alphaSum=0;
		for (int i=0; i<m; i++)
			out+=mylgamma(alpha[i]);
		for (int i=0; i<m; i++)
			alphaSum+=alpha[i];
		out-=mylgamma(alphaSum);
		if (logged) return out;
		else return exp(out);
	}
	static double lbeta (std::vector<double> & alpha, bool logged=true) {
        int m=alpha.size();   
		double out=0;
		double alphaSum=0;
		for (int i=0; i<m; i++)
			out+=mylgamma(alpha[i]);
		for (int i=0; i<m; i++)
			alphaSum+=alpha[i];
		out-=mylgamma(alphaSum);
		if (logged) return out;
		else return exp(out);
	}
	static double lbeta (int* alpha, int m, bool logged=true) {
		double out=0;
		int alphaSum=0;
		for (int i=0; i<m; i++)
			out+=mylgamma(alpha[i]);
		for (int i=0; i<m; i++)
			alphaSum+=alpha[i];
		out-=mylgamma(alphaSum);
		if (logged) return out;
		else return exp(out);
	}
	
	// given the sequence no. of a cell from the lower left triangle, return the row and column
	// row will always be greater or equal to col
	// seqno, row, and col are all 0-based
	static void getRowCol (int seqno, int & row, int & col) {
		seqno ++; // make it start from 1
		row = floor (sqrt ((double)(seqno * 2)) );
		col = seqno - row * (row + 1) /2 - 1; 
		if (col >= 0) 
			row ++;
		else 
			col += row;
//		// code for testing
//		int counter=-1;
//		for (int j1=1; j1<10; j1++)
//			for (int j2=0; j2<j1; j2++) {
//				counter ++;
//				int row, col;
//				MathPlus::getRowCol(counter, row, col);
//				cout << j1-row << ",    " << j2-col << endl;
//			}    
	}

};
#endif
