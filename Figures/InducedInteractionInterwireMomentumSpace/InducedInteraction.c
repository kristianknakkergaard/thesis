#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
/*Units: */
/*Energy: the Fermi energy \epsilon_F */
/*Momentum: the fermi momentum k_F */
/*Length: */
/*All is for T=0*/

/*nB/nF^3, kF * aB HERE*/
/*SET mB/mF, nB/nF^3 and kF * aBF */
double VFF12(double q, double d, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 

  double factor = 16.0/(M_PI * M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);
  double f      = - factor * gsl_sf_bessel_K0(sqrt(pow(q * d, 2.0) + pow(2.0 * d/xi, 2.0)  ));
  return f;
}

int
main (void)
{
  /*variables */
  
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double rBB = 0.01;
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB = 100.0;
  double d, d_low = 0.6, d_up = 2.0, dd = 0.4;

  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);
  double retneg = 4.0/(M_PI * M_PI) * pow(1.0/mB, 2.0) * nB * aB;

  /*q's*/
  double q_low = -10.0, q_up = 10.0, dq = 0.01; 
  
  printf("%s \t %s \t %s \n", "q", "Induced Interaction", "retneg" );
  printf("\n\n");

  for (d = d_low; d < d_up; d+= dd)
  {
    for ( double q = q_low; q < q_up; q += dq )
    {
      if (retneg < 1.0)
      {
        break;
      }
      printf("%lg \t %lg \t %lg \n", q, VFF12(q, d, rBB, rBF, nB, mB), retneg );
    }
    printf("\n\n");
  }


  return 0;
}