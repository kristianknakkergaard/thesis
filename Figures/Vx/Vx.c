#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

double Vx (double x, double lt, double xi)
{
  double Fx = lt/xi + fabs(x)/(sqrt(2.0) * lt);
  return - sqrt(M_PI/2.0) * 1.0/lt * exp(-sqrt(2)*fabs(x)/xi) * exp(Fx * Fx) * erfc(Fx); 
}

double Vx_asymp (double x, double xi)
{
  return - exp(-sqrt(2)*fabs(x)/xi) / fabs(x); 
}

int
main (void)
{
  double BBscatlength = 0.3;
  double densityratio = 1;
  double xi = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/
  double lt, lt_low = 0.03, lt_high = 0.2, dlt = 0.03;
  double x, x_low = -1, x_high = 1, dx = 0.001; 

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
    for (x = x_low; x<= x_high; x+=dx)
    {
      printf("%lg \t %lg \t %lg \t %lg \n", lt, x, Vx(x,lt,xi), Vx_asymp(x,xi));
    }
    printf("\n\n");
  }
  
  return 0;
}