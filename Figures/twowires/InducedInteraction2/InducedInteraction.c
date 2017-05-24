#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
/*Units: */
/*Energy: the Fermi energy \epsilon_F */
/*Momentum: the fermi momentum k_F */
/*Length: */
/*All is for T=0*/

/*nB/nF^3 HERE*/
/*SET mB/mF, nB/nF^3 and kF * aBF */
double WFF11(double k, double kprime, double rBB, double rBF, double nB, double mB)
{
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/ 

  double factor = 4.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;
  double f      = - factor * log( ( pow( k + kprime, 2.0 ) + 2.0 / (xi * xi) ) / ( pow( k - kprime, 2.0 ) + 2.0 / (xi * xi) ) );
  return f; 
}

double VFF12(double q, double d, double rBB, double rBF, double nB, double mB)
{
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/

  double factor = 16.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;
  double f      = - factor * gsl_sf_bessel_K0( sqrt( pow(q * d, 2.0) + 2.0 * pow( d / xi, 2.0) ) );
  return f;
}

double WFF12(double k, double kprime, double d, double rBB, double rBF, double nB, double mB)
{
  return 1.0 / 2.0 * ( VFF12(k + kprime, d, rBB, rBF, nB, mB) + VFF12(k - kprime, d, rBB, rBF, nB, mB) ); 
}

int
main (void)
{
  /*variables */
  double rBF  = 0.11;    /*(nB * aBF^3)^(1/3) */
  double rBB  = 0.01;
  double nB   = 100.0;
  double mB   = 7.0/40.0; /*mB/mF*/
  
  /*distances:*/
  const int Nd        = 5;
  double d_string[Nd] = {0.720, 0.735, 0.750, 0.765, 0.775};
  double d;

  /*k-values:*/
  double k_low = -10.0, k_up = 10.0, dk = 0.01;

  double W11, W12;

  printf("%s \t %s \t %s \t %s \n", "k", "Wind11", "Wind12", "d");
  
  for (int id = 0; id < Nd; ++id)
  {
    d = d_string[id];
    for (double k = k_low; k <= k_up; k += dk )
    {
      W11 = WFF11(k, 1.8,    rBB, rBF, nB, mB);
      W12 = WFF12(k, 0.0, d, rBB, rBF, nB, mB);

      printf("%lg \t %lg \t %lg \t %lg \n", k, W11, W12, d);
    }
    printf("\n \n");
  }


  return 0;
}