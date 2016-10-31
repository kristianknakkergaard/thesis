#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>

/*units: 
momentum: k_F 
distance: 1/k_F*/

double Vq (double q, double lt, double xi)
{
  double F = lt*lt/2.0 * (q * q + 2/(xi * xi));
  return -gsl_sf_exp(F) * gsl_sf_expint_E1(F); 
}

int
main (void)
{
  double BBscatlength = 0.3;
  double densityratio = 1;
  double xi = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/
  double lt, lt_low = 0.03, lt_high = 0.2, dlt = 0.03;
  double q, q_low = -20, q_high = 20, dq = 0.001; 

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
    for (q = q_low; q<= q_high; q+=dq)
    {
      printf("%lg \t %lg \t %lg \n", lt, q, Vq(q,lt,xi));
    }
    printf("\n\n");
  }
  
  return 0;
}