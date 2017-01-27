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
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/ 
  double factor = 8.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;

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