// Author: Krisztian Sebestyen ksebestyen@gmail.com

//#include "matrix.h"
#include <float.h> //DBL_EPSILON
#include <R_ext/Lapack.h>
#include <R.h>
//#include <R_ext/Applic.h>
//#include <R_ext/BLAS.h>
//#include <R_ext/RS.h> //R definitions for 'extending' R, registering functions,...

#define PRINTF Rprintf
#define	MAX(A,B)	((A) > (B) ? (A) : (B))
#define	MIN(A,B)	((A) < (B) ? (A) : (B))


// LEADING DIMENSION 'LDA' or 'LD', comments collected from the 'net'
// 'ld' : 'leading dimension of matrix' is equal to the number of elements in its major dimension.
//      : the distance between two neighboring elements along minor dimension.
// It is used to access a sub-matrix of a larger matrix
// The leading dimension for a two-dimensional array is an increment that is used to find the 
// starting point for the matrix elements in each successive column of the array. 
// The leading dimension must always be positive. It must always be greater than or equal to the
// minor dimension of the matrix.
#define CX(i,j,ld) (((j)*(ld))+(i)) //indexing matrices in column-major order
#define RX(i,j,ld) (((i)*(ld))+(j)) //indexing matrices in row-major order
typedef enum {COLUMN_MAJOR = 0,ROW_MAJOR = 1}MAJOR_ORDER;



// matrix multiplication D %*% X %*% D2 for D,D2 diagonal matrices and X square matrix in column-major order
// x,y are square matrices in column order, d1 and d2 represent the diagonals of diagonal matrices
void dxd_(int* _n, double* d1,double* x,double* d2,double* result){
	int i,j;
	int n = *_n;
	for(j = 0;j < n;j++)
		for(i = 0;i < n;i++)
			result[j*n + i] = d1[i] * x[j*n + i] * d2[j];
}

void dxd(int* n, double* d1,double* _x,double* d2,double* _y){
	int i,j;
	double* x = _x;
	double* y = _y;
	for(j = 0;j < *n;j++){
		for(i = 0;i < *n;i++){
			*y = d1[i] * (*x) * d2[j];
			y++;x++;
		}
	}
}


void symprod(int* nrow,int* ncol,double* S,double* X,double* Y)
{

	char*  	SIDE = "L";
	char*  	UPLO = "U";
	double alpha = 1.0;
	double beta = 0.0;
	
	F77_CALL(dsymm)( 	
		SIDE,UPLO,
		nrow,ncol,
		&alpha,
		S,nrow,
		X,nrow,
		&beta,
		Y,nrow 
	); 
	
}


void txSy(int* n,double* S,double* x,double* y,double* temp,double* out){

	double alpha = 1.0;
	double beta = 0.0;
	int ione = 1;
	F77_CALL(dsymv)("U",n,&alpha,S,n,y,&ione,&beta,temp,&ione); // temp = S %*% y
	*out = F77_CALL(ddot)(n, x, &ione, temp, &ione); // out = x %*% temp
}


