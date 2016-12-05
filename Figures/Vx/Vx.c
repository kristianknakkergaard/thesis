#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

double Vx (double x, double lt, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 
  double factor = 8.0/(M_PI * M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);

  double Fx = lt/xi + fabs(x)/(sqrt(2.0) * lt);
  return - factor * sqrt(M_PI/2.0) * 1.0/lt * exp(-sqrt(2)*fabs(x)/xi) * exp(Fx * Fx) * erfc(Fx); 
}

double Vx_asymp (double x, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 
  double factor = 8.0/(M_PI * M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);

  return - factor * exp(-sqrt(2)*fabs(x)/xi) / fabs(x); 
}

int
main (void)
{
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double rBB = 0.01;
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB = 100.0;

  double lt, lt_low = 0.05, lt_high = 0.2, dlt = 0.05;
  double x, x_low = -5, x_high = 5, dx = 0.001; 

  for (lt = lt_low; lt <= lt_high; lt+=dlt)
  {
    for (x = x_low; x<= x_high; x+=dx)
    {
      printf("%lg \t %lg \t %lg \t %lg \n", lt, x, Vx(x, lt, rBB, rBF, nB, mB), Vx_asymp(x, rBB, rBF, nB, mB));
    }
    printf("\n\n");
  }
  
  return 0;
}