#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
/*Units: */
/*Energy: the Fermi energy \epsilon_F */
/*Momentum: the fermi momentum k_F */
/*Length: */
/*All is for T=0*/

/*nB/nF^3, kF * aB HERE*/
/*SET mB/mF, nB/nF^3 and kF * aBF */
double WFF0(double k, double kprime, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 

  double factor = 4.0/(M_PI * M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);
  double f      = - factor * log( ( pow( k + kprime, 2.0) + 2.0/(xi * xi) )/( pow( k - kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}

int
main (void)
{
  /*variables */
  
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double rBB = 0.01;
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB;

  /*nB:*/
  double Bdist, Bdist_low = 0.0001, Bdist_up = 0.8, dBdist = 0.0001;

  double aB, retneg;

  
  printf("%s \t %s \t %s \n", "Bdist", "Induced Interaction", "retneg" );
  printf("\n\n");

  for ( Bdist = Bdist_low; Bdist < Bdist_up; Bdist += dBdist )
  {
    nB     = pow(Bdist, - 3.0);
    aB     = M_PI * rBB / pow(nB, 1.0/3.0);
    retneg = 4/(M_PI * M_PI) * pow(1.0/mB, 2.0) * nB * aB;

    /*if (retneg < 5.0)
    {
      break;
    }*/
    
    printf("%lg \t %lg \t %lg \n", Bdist, WFF0(1.0, 1.0, rBB, rBF, nB, mB), retneg );
  }


  return 0;
}