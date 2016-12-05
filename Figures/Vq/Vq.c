#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>

/*units: 
momentum: k_F 
distance: 1/k_F*/

double Vq (double q, double lt, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 
  double factor = 8.0/(M_PI * M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);

  double F = lt*lt/2.0 * (q * q + 2/(xi * xi));
  return - factor * gsl_sf_exp(F) * gsl_sf_expint_E1(F); 
}

int
main (void)
{
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double rBB = 0.01;
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB = 100.0;

  double lt, lt_low = 0.05, lt_high = 0.2, dlt = 0.05;
  double q, q_low = -20, q_high = 20, dq = 0.001; 

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
    for (q = q_low; q<= q_high; q+=dq)
    {
      printf("%lg \t %lg \t %lg \n", lt, q, Vq(q,lt, rBB, rBF, nB, mB));
    }
    printf("\n\n");
  }
  
  return 0;
}