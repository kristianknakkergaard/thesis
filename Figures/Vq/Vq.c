#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>


double Vq (double q, double lt)
{
  double F = lt*lt/2.0 * (q*q + 2);
  return -gsl_sf_exp(F) * gsl_sf_expint_E1(F); 
}

int
main (void)
{
  double lt, lt_low = 0.1, lt_high = 2, dlt = 0.3;
  double q, q_low = -10, q_high = 10, dq = 0.001; 

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
    for (q = q_low; q<= q_high; q+=dq)
    {
      printf("%lg \t %lg \t %lg \n", lt, q, Vq(q,lt));
    }
    printf("\n\n");
  }
  
  return 0;
}