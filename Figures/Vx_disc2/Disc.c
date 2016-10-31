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

double Disc (double x, double lt, double xi, double disc)
{
  double Fx = lt/xi + x/(sqrt(2.0) * lt);
  return sqrt(M_PI/2.0) * x/lt * erfcx(Fx) - disc; 
}


int
main (void)
{
  double disc = 0.95;
  double BBscatlength = 0.3;
  double densityratio = 1;
  double xi = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/

  double x_low=0.001, x_up=300, dx = 0.001;
  double lt_0 = 0.001, lt_up = 0.3, dlt = 0.001;

  for (double lt = lt_0; lt < lt_up; lt+=dlt)
  {
    for (double x = x_low; x < x_up; x+=dx)
    {
      if (Disc(x,lt,xi,disc) >= 0 )
      {
        printf("%lg \t %lg \t %lg \t %lg \n", lt, x, Disc(x,lt,xi,disc), 1.0/sqrt(1-disc)*lt );
        break;
      }
    }
  }


  return 0;
}