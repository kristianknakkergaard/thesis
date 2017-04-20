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

double Vindx (double x, double rBB, double rBF, double nB, double mB, double eps0, double r)
{
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/
  double factor = 8.0 / M_PI * (mB + 1.0 / mB + 2.0) * eps0 * pow(nB, 1.0 / 3.0) * rBF * rBF;

  return - factor * exp(- sqrt(2.0) * fabs(x) / xi ) / pow( fabs(x), r); 
}

double WFF0(double k, double kprime, double rBB, double rBF, double nB, double mB, int N, double eps0, double r)
{
  double W = 0.0;

  for (int i = 1; i <= N; ++i)
  {
    double x  = (double) i;
    W += 1.0 / 2.0 * ( cos((kprime - k) * i ) - cos( (kprime + k) * i ) ) * Vindx( x, rBB, rBF, nB, mB, eps0, r); 
  }

  return W; 
}

double Deltaguess(double k)
{
  return - 0.2 * sin(2.0 * k);
}


int
main (void)
{
  /*variables */
  double rBB     = 0.01; /*(nB * aBB^3)^(1/3) < 0.03 atmost! (1 percent depletion) */
  double rBF     = 0.15;    /*(nB * aBF^3)^(1/3) */
  double nB      = 100.0; 
  double mB      = 7.0/40.0; /*mB/mF*/
  
  double t1      = 1.0;  /*the unit of energy, also gives unit of length */
  double t2      = 1.0;
  double eps0    = 3.0; 
  
  int    N       = 400; 
  double Ndouble = (double)N;
  double r       = 1.0;

  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  
  double gapintegrand;
  double gapintegral;
  double D_maxk;

  double k;
  double kprime;
  int check = 0;


  double xi       = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0); 
  double vFoverc0 = sqrt(M_PI) / 2.0 * mB * 1.0 / sqrt(rBB) * pow(nB, -1.0 / 3.0); 
  double strength = 8.0 / M_PI * (mB + 1.0 / mB + 2.0) * eps0 * pow(nB, 1.0 / 3.0) * rBF * rBF;

  fprintf(stderr, "strength = %lg, (nBaB^3)^(1/3) = %lg, (nBaBF^3)^(1/3) = %lg, mB/mF = %lg, nB/nF^3 = %lg, kF*xi = %lg, vF/c0 = %lg \n", strength, rBB, rBF, mB, nB, xi, vFoverc0);
  fprintf(stderr, "\n \n");

  /*vectors and matrices:*/
  gsl_vector *D          = gsl_vector_calloc(N);
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);


  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set( D, i, Deltaguess(k) );

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, rBB, rBF, nB, mB, N, eps0, r) );
    }
  }

  double mu = 0.5;
  double muintegral = 0.0, muintegrand;
  
  
  /*T = 0: */
  for (int j = 0; j < 1000; ++j)
  {
    D_maxk = gsl_vector_max(D);
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gapintegral = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
        Deltakprime   = gsl_vector_get(D, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
        gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
        gapintegral  += gapintegrand; 
      }
      gsl_vector_set(D, i, gapintegral);
    }
  
    if ( fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-4 )
    {
      break;
    }

    if (gsl_vector_max(D) < 1e-4 )
    {
      break;
    }
    ++check; 
  }

  /*variables for CS1 */
  double CS1 = 0;
  double CS1integrand, dDelta;
  double depsilon; 

  for (int iprime = 0; iprime < N; ++iprime)
  {
    kprime        = ((double)iprime) * dk + k_low;
    epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
    Deltakprime   = gsl_vector_get(D, iprime);
    EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
    muintegrand   = 1.0 / (2.0 * Ndouble) * (1.0 -  epsilonkprime / EFkprime);
    muintegral   += muintegrand;
    
    if (iprime < N - 1)
      {
        dDelta       = (gsl_vector_get(D, iprime + 1) - Deltakprime ) / dk;
        depsilon     = t1 * sin(kprime) + 2.0 * t2 * sin(2.0 * kprime);
        CS1integrand = 1.0 / (2.0 * M_PI) * (Deltakprime * depsilon - epsilonkprime * dDelta)/(EFkprime * EFkprime);
        CS1         += dk * CS1integrand; 
      }
  }

  /*We list the function values:*/
  printf("%s \t %s \t %s \t %s \t %s \t %s \n", "k", "Deltak", "mu", "Nmu/N", "check", "CS1" );
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \t %lg \t %i \t %lg \n", k, gsl_vector_get(D, i), mu, muintegral, check, CS1);
  }
  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  return 0;
}