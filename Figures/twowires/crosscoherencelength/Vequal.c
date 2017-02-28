#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

double Vequal (double d, double rBB, double nB, double mB, double x0, double ratio)
{
  double xi    = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/

  double alpha = sqrt(1.0 + pow(d / x0, 2.0));
  double f = alpha * exp(-sqrt(2.0) * x0 / xi * (1.0 - alpha ));
  return f - ratio;
}


int
main (void)
{
  double rBB   = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double mB    = 7.0/40.0; /*mB/mF*/
  double ratio = 2.0;
  double x0    = 1.0/sqrt(3.0);
  double nB    = 0.0;

  double Bdist_low = 0.0001, Bdist_up = 10.0, dBdist = 0.001;
  double d_low  = 0.00001, d_up  = 1.0, dd = 0.0001;

  int Bdist_iter = 0; 
  for (double Bdist = Bdist_low; Bdist < Bdist_up; Bdist+=dBdist)
  {
    nB = pow(Bdist, -3.0);
    for (double d = d_low; d < d_up; d+=dd)
    {
      if ( Vequal (d, rBB, nB, mB, x0, ratio) >= 0.0 )
      {
        printf("%lg \t %lg \t %lg \n", Bdist, d, Vequal (d, rBB, nB, mB, x0, ratio));
        break;
      }
    }
    fprintf(stderr, "Bdist_iter = %i\n", ++Bdist_iter);
  }

  return 0;
}