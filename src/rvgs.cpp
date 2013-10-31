// modified by Youyi Fong

/* -------------------------------------------------------------------------- 
 * This is an ANSI C library for generating random variates from six discrete 
 * distributions
 *
 *      Generator         Range (x)     Mean         Variance
 *
 *      Bernoulli(p)      x = 0,1       p            p*(1-p)
 *      Binomial(n, p)    x = 0,...,n   n*p          n*p*(1-p)
 *      Equilikely(a, b)  x = a,...,b   (a+b)/2      ((b-a+1)*(b-a+1)-1)/12
 *      Geometric(p)      x = 0,...     p/(1-p)      p/((1-p)*(1-p))
 *      Pascal(n, p)      x = 0,...     n*p/(1-p)    n*p/((1-p)*(1-p))
 *      Poisson(m)        x = 0,...     m            m
 * 
 * and seven continuous distributions
 *
 *      Uniform(a, b)     a < x < b     (a + b)/2    (b - a)*(b - a)/12 
 *      Exponential(m)    x > 0         m            m*m
 *      Erlang(n, b)      x > 0         n*b          n*b*b
 *      Normal(m, s)      all x         m            s*s
 *      Lognormal(a, b)   x > 0            see below
 *      Chisquare(n)      x > 0         n            2*n 
 *      Student(n)        all x         0  (n > 1)   n/(n - 2)   (n > 2)
 *
 * For the a Lognormal(a, b) random variable, the mean and variance are
 *
 *                        mean = exp(a + 0.5*b*b)
 *                    variance = (exp(b*b) - 1) * exp(2*a + b*b)
 *
 * Name              : rvgs.c  (Random Variate GeneratorS)
 * Author            : Steve Park & Dave Geyer
 * Language          : ANSI C
 * Latest Revision   : 10-28-98
 * --------------------------------------------------------------------------
 */

#include "math.h"

#include "rvgs.h"

using namespace std;

const double EE = 2.718281828459046;

   long rBernoulli(double p)
/* ========================================================
 * Returns 1 with probability p or 0 with probability 1 - p. 
 * NOTE: use 0.0 < p < 1.0                                   
 * ========================================================
 */ 
{
  return ((myunif_rand() < (1.0 - p)) ? 0 : 1);
}

   long rBinomial(long n, double p)
/* ================================================================ 
 * Returns a binomial distributed integer between 0 and n inclusive. 
 * NOTE: use n > 0 and 0.0 < p < 1.0
 * ================================================================
 */
{ 
  long i, x = 0;

  for (i = 0; i < n; i++)
    x += rBernoulli(p);
  return (x);
}

long rUniform(long a, long b)
/* ===================================================================
 * Returns an equilikely distributed integer between a and b inclusive. 
 * NOTE: use a < b
 * ===================================================================
 */
{
  return (a + (long) ((b - a + 1) * myunif_rand()));
}
int rUniform(int a, int b) {
  return (a + (int) ((b - a + 1) * myunif_rand()));
}

   long rGeometric(double p)
/* ====================================================
 * Returns a geometric distributed non-negative integer.
 * NOTE: use 0.0 < p < 1.0
 * ====================================================
 */
{
  return ((long) (log(1.0 - myunif_rand()) / log(p)));
}

   long rPascal(long n, double p)
/* ================================================= 
 * Returns a Pascal distributed non-negative integer. 
 * NOTE: use n > 0 and 0.0 < p < 1.0
 * =================================================
 */
{ 
  long i, x = 0;

  for (i = 0; i < n; i++)
    x += rGeometric(p);
  return (x);
}

   long rPoisson(double m)
/* ================================================== 
 * Returns a Poisson distributed non-negative integer. 
 * NOTE: use m > 0
 * ==================================================
 */
{ 
  double t = 0.0;
  long   x = 0;

  while (t < m) {
    t += rExponential(1.0);
    x++;
  }
  return (x - 1);
}

   double rUniform(double a, double b)
/* =========================================================== 
 * Returns a uniformly distributed real number between a and b. 
 * NOTE: use a < b
 * ===========================================================
 */
{ 
  return (a + (b - a) * myunif_rand());
}

   double rExponential(double m)
/* =========================================================
 * Returns an exponentially distributed positive real number. 
 * NOTE: use m > 0.0
 * =========================================================
 */
{
  return (-m * log(1.0 - myunif_rand()));
}

   double rErlang(long n, double b)
/* ================================================== 
 * Returns an Erlang distributed positive real number.
 * NOTE: use n > 0 and b > 0.0
 * ==================================================
 */
{ 
  long   i;
  double x = 0.0;

  for (i = 0; i < n; i++) 
    x += rExponential(b);
  return (x);
}

   double rNormal(double m, double s)
/* ========================================================================
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Uses a very accurate approximation of the normal idf due to Odeh & Evans, 
 * J. Applied Statistics, 1974, vol 23, pp 96-97.
 * ========================================================================
 */
{ 
  const double p0 = 0.322232431088;     const double q0 = 0.099348462606;
  const double p1 = 1.0;                const double q1 = 0.588581570495;
  const double p2 = 0.342242088547;     const double q2 = 0.531103462366;
  const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850;
  const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2;
  double u, t, p, q, z;

  u   = myunif_rand();
  if (u < 0.5)
    t = sqrt(-2.0 * log(u));
  else
    t = sqrt(-2.0 * log(1.0 - u));
  p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
  q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
  if (u < 0.5)
    z = (p / q) - t;
  else
    z = t - (p / q);
  return (m + s * z);
}

   double rLognormal(double a, double b)
/* ==================================================== 
 * Returns a lognormal distributed positive real number. 
 * NOTE: use b > 0.0
 * ====================================================
 */
{
  return (exp(a + b * rNormal(0.0, 1.0)));
}

   double rChisquare(long n)
/* =====================================================
 * Returns a chi-square distributed positive real number. 
 * NOTE: use n > 0
 * =====================================================
 */
{ 
  long   i;
  double z, x = 0.0;

  for (i = 0; i < n; i++) {
    z  = rNormal(0.0, 1.0);
    x += z * z;
  }
  return (x);
}

   double rStudent(long n)
/* =========================================== 
 * Returns a student-t distributed real number.
 * NOTE: use n > 0
 * ===========================================
 */
{
  return (rNormal(0.0, 1.0) / sqrt(rChisquare(n) / n));
}

unsigned long rMultnomial(double *adPr, unsigned long k)
//return one multinomial rv bewteen 0 and k-1
{
    double dP = myunif_rand();
    double dSum = 0.0; 
    unsigned long i = 0;

    //assert(adPr!=NULL);

    for(i=0; i<k; i++)
    {
        dSum = dSum + adPr[i];
        if(dSum > dP)
        {
            return i;
        }
    }

    return k-1;
}
unsigned short rMultnomial(double *adPr, unsigned short k)
//return one multinomial rv bewteen 0 and k-1
{
    double dP = myunif_rand();
    double dSum = 0.0; 
    unsigned short i = 0;

    //assert(adPr!=NULL);

    for(i=0; i<k; i++)
    {
        dSum = dSum + adPr[i];
        if(dSum > dP)
        {
            return i;
        }
    }

    return k-1;
}

double rGamma(double dAl, bool give_log)
{
    double dR1,dR2,dAA,dX,dW,dC1,dC2,dC3,dC4,dC5;
    double logdX;
    if(dAl<=0.0)
    {
        return log(0.0);
    }
    if(dAl == 1.0)
    {
		double r=rExponential(1.0);
		return give_log?log(r):r;
    }
    if(dAl<1.0)
    {
        dAA=(dAl+EE)/EE;
        do
        {
            dR1=myunif_rand();
            dR2=myunif_rand();
            if(dR1>1.0/dAA)
            {
                dX = -log(dAA*(1.0-dR1)/dAl);
                if(dR2<pow(dX,(dAl-1.0)))
                {
                    return give_log?log(dX):dX;
                }
            }
            else
            {
                logdX = log(dAA*dR1)*(1.0/dAl);
                dX=exp(logdX);
                //dX = pow((dAA*dR1),(1.0/dAl));
                if (dR2<exp(-dX))
                {
                    return give_log?logdX:dX;
                }
            }
        }while(dR2<2);
    }
    else
    {
        dC1=dAl-1.0;
        dC2=(dAl-1.0/(6.0*dAl))/dC1;
        dC3=2.0/dC1;
        dC4=dC3+2.0;
        dC5=1.0/sqrt(dAl);
        do
        {
            do
            {
                dR1=myunif_rand();
                dR2=myunif_rand();
                if(dAl>2.5)
                {
                    dR1=dR2+dC5*(1.0-1.86*dR1);
                }
            }while(dR1<=0 || dR1 >= 1);
            dW=dC2*dR2/dR1;
            if(dC3*dR1+dW+1/dW <= dC4)
            {
                return give_log?log(dC1*dW):dC1*dW;
            }
            if(dC3*log(dR1)-log(dW)+dW<1)
            {
                return give_log?log(dC1*dW):dC1*dW;
            }
        }while(dR2<2);
    }
	return log(0.0);
}
