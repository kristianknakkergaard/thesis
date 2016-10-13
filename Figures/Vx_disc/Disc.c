#include <stdlib.h>
#include <math.h>
#include <tgmath.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

double erfcx(double x)
{
  if (x < 25)
  {
    return exp(x * x) * erfc(x) ;
  }
  else
  {
    double y = 1.0 / x;
    double z = y * y;
    double s = y * (1.0 + z * (-0.5 + z * (0.75 + z * (-1.875 + z * (6.5625 - 29.53125 * z)))));
    return s * 0.564189583547756287;
  }
}

double Disc (double x, double lt, double disc)
{
  double Fx = lt + x/(sqrt(2.0) * lt);
  return sqrt(M_PI/2.0) * x/lt * erfcx(Fx) - disc; 
}


int
main (void)
{
  double disc = 0.95;

  double x_low, x_up, dx = 0.001;
  double lt_0 = 0.001, lt_up = 2, dlt = 0.01;

  for (double lt = lt_0; lt < lt_up; lt+=dlt)
  {
    if (lt < 0.25)
    {
      x_low = 0.001;
      x_up = 2.25;
    }
    else 
    {
      x_low = 20*lt*lt;
      x_up = 40*lt*lt;
    }
    for (double x = x_low; x < x_up; x+=dx)
    {
      if (Disc(x,lt,disc) >= 0 )
      {
        printf("%lg \t %lg \t %lg \t %lg \n", lt, x, Disc(x,lt,disc), 27.0*lt*lt );
        break;
      }
    }
  }


  return 0;
}