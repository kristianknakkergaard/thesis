#include <stdlib.h>
#include <tgmath.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_exp.h>


double Vq0 (double lt)
{
  double F = lt*lt;
  return -gsl_sf_exp(F) * gsl_sf_expint_E1(F); 
}

int
main (void)
{
  double lt, lt_low = 0.0001, lt_high = 1, dlt = 0.00001;

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
      printf("%lg \t \t %lg \t \t %lg \n", 1.0/lt, Vq0(lt), 1.93*log(lt));
  }
  
  return 0;
}