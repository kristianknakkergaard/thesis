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
  
  double rBF = 0.11;    /*(nB * aBF^3)^(1/3) */
  double rBB = 0.01;
  double d   = 0.7;
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB;

  /*nB:*/
  double Bdist_low = 0.0005, Bdist_up = 1.5, dBdist = 0.0005;

  double vFoverc0;

  
  printf("%s \t %s \t %s \t %s \n", "Bdist", "Wind11", "Wind12", "vFoverc0" );
  
  for ( double Bdist = Bdist_low; Bdist < Bdist_up; Bdist += dBdist )
  {
    nB       = pow(Bdist, - 3.0);
    vFoverc0 = sqrt(M_PI) / 2.0 * mB * 1.0 / sqrt(rBB) * pow(nB, -1.0 / 3.0);
    
    printf("%lg \t %lg \t %lg \t %lg \n", Bdist, WFF11(1.0, 1.0, rBB, rBF, nB, mB), WFF12(1.0, 1.0, d, rBB, rBF, nB, mB), vFoverc0 );
  }


  return 0;
}