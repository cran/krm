// modified by Youyi Fong

/* ------------------------------------------------------------- 
 * Name            : rvgs.h (header file for the library rvgs.c)
 * Author          : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 11-03-96
 * -------------------------------------------------------------- 
 */

#if !defined( _RVGS_ )
#define _RVGS_


#define AA_KINDS 20

extern "C" double myunif_rand();

long rBernoulli(double p);
long rBinomial(long n, double p);
long rUniform(long a, long b);
int rUniform(int a, int b);
long rGeometric(double p);
long rPascal(long n, double p);
long rPoisson(double m);
// added by Youyi Fong
unsigned short rMultnomial(double *adPr, unsigned short k);
unsigned long  rMultnomial(double *adPr, unsigned long k);

double rUniform(double a, double b);
double rExponential(double m);
double rErlang(long n, double b);
double rNormal(double m, double s);
double rLognormal(double a, double b);
double rChisquare(long n);
double rStudent(long n);
// added by Youyi Fong
// the following Gamma rv has rate 1, for rate beta, multiply by 1/beta
double rGamma(double dAl, bool give_log=false);


#endif

