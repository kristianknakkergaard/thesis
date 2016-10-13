#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

double Vx (double x, double lt)
{
  double Fx = lt + fabs(x)/(sqrt(2.0) * lt);
  return - sqrt(M_PI/2.0) * 1.0/lt * exp(-sqrt(2)*fabs(x)) * exp(Fx * Fx) * erfc(Fx); 
}

double Vx_asymp (double x)
{
  return - exp(-sqrt(2)*fabs(x)) / fabs(x); 
}

int
main (void)
{
  double lt, lt_low = 0.1, lt_high = 2, dlt = 0.3;
  double x, x_low = -10, x_high = 10, dx = 0.001; 

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
    for (x = x_low; x<= x_high; x+=dx)
    {
      printf("%lg \t %lg \t %lg \t %lg \n", lt, x, Vx(x,lt), Vx_asymp(x));
    }
    printf("\n\n");
  }
  
  return 0;
}