#include <stdlib.h>
#include <math.h>
#include <tgmath.h>

double Vx (double x, double lt, double rBB, double rBF, double nB, double mB)
{
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/ 
  double factor = 8.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;
  double Fx     = lt / xi + fabs(x)/(sqrt(2.0) * lt);

  return - factor * sqrt(M_PI / 2.0) * 1.0 / lt * exp( -sqrt(2) * fabs(x) / xi ) * exp( Fx * Fx ) * erfc( Fx ); 
}

double Vx_asymp (double x, double rBB, double rBF, double nB, double mB)
{
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/
  double factor = 8.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;

  return - factor * exp(-sqrt(2)*fabs(x)/xi) / fabs(x); 
}

int
main (void)
{
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double rBB = 0.01;
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB = 100.0;

  double lt, lt_low = 0.2, lt_high = 0.7, dlt = 0.2;
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