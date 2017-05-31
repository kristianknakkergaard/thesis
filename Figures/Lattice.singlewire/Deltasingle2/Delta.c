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

double Vindx (double x, double xi, double strength)
{
  return - strength * exp( - sqrt(2.0) * fabs(x) / xi );
}

double WFF0( double k, double kprime, double xi, double strength, int N )
{
  double W = 0.0;

  for (int i = 1; i < N; ++i)
  {
    double x  = (double) i;
    W += 1.0 / 2.0 * ( cos((kprime - k) * i ) - cos( (kprime + k) * i ) ) * Vindx( x, xi, strength ); 
  }

  return W; 
}

double Deltaguess(double k)
{
  return 0.0 * sin(k) + 1.0 * sin(2.0 * k);
}


int
main (void)
{
  /*restricted variables */
  double t1       = 1.0;  
  double strength = 2.0; /*total strength of interaction*/
  double lt       = 0.05;

  /*free parameters:*/
  double xi       = 5.0; /*we will eventually want to control the range rather than the Bose gas parameter*/
  double t2       = 1.0;
  int    N        = 200; 
  double Ndouble  = (double)N;
  double mu       = 0.0;
  double filling;

  /*derived parameters:*/
  double rBB      = pow(lt / (2.0 * xi), 2.0); /*(nB * aBB^3)^(1/3) < 0.03 atmost! (1 percent depletion) */

  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  
  double gapintegrand;
  double gapintegral;
  double muintegral, muintegrand;
  double D_maxk;
  double Deltai;

  double k;
  double kprime;
  int check = 0;

  /*iterations and convergence*/
  double convergence = 1e-7;
  int    iterations  = 1000;

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
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, xi, strength, N) );
    }
  }

  /*T = 0: */
  for (int j = 0; j < iterations; ++j)
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

    muintegral  = 0.0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
      Deltakprime   = gsl_vector_get( D, iprime );
      EFkprime      = pow( epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0 );
      muintegrand   = 1.0 / ( 2.0 * Ndouble ) * ( 1.0 -  epsilonkprime / EFkprime );
      muintegral   += muintegrand;
    }

    if ( fabs(gsl_vector_max(D) - D_maxk)/D_maxk < convergence && fabs(muintegral - filling) < convergence)
    {
      break;
    }
    if ( fabs(gsl_vector_max(D)) < convergence && fabs(muintegral - filling) < convergence )
    {
      break;
    }

    filling = muintegral;
    fprintf(stderr, "Deltamax = %lg, \t n = %lg,  \t check = %i \n", gsl_vector_max(D), filling, ++check );
  }

  /*variables for CS1 */
  double CS1 = 0;
  double kprimenext;
  double epsilonkprimenext; 
  int Nhalf = (int) Ndouble / 2.0; 
  

  for (int iprime = 0; iprime < Nhalf; ++iprime)
  {
    kprime            = ((double)iprime      ) * dk + k_low;
    kprimenext        = ((double)iprime + 1.0) * dk + k_low;

    epsilonkprime     = - t1 * cos(kprime)     - t2 * cos(2.0 * kprime)     - mu;
    epsilonkprimenext = - t1 * cos(kprimenext) - t2 * cos(2.0 * kprimenext) - mu;

    if (epsilonkprime < 0 && epsilonkprimenext > 0)
    {
      CS1 += copysign(1.0, gsl_vector_get(D, iprime)); 
    }
    if (epsilonkprime > 0 && epsilonkprimenext < 0)
    {
      CS1 -= copysign(1.0, gsl_vector_get(D, iprime)); 
    }
  }

  /*Ground state energy: */

  double E0integrand; 
  double E0integral = 0; 
  for (int iprime = 0; iprime < N; ++iprime)
  {
    kprime        = ((double)iprime) * dk + k_low;
    epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
    Deltakprime   = gsl_vector_get(D, iprime);
    EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

    E0integrand   = - 1.0 / (4.0 * Ndouble) * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
    E0integral   += E0integrand; 
  }

  fprintf(stderr, "E0 + filling * mu = %lg, mu = %lg, filling = %lg, nBaB = %lg, xi / a = %lg, strength = %lg, CS1 = %lg \n", E0integral + filling * mu, mu, filling, rBB, xi, strength, CS1);
  fprintf(stderr, "\n \n");


  /*We list the function values:*/ 
  printf("%s \t %s \t %s \t %s \t %s \t %s \n", "k", "Deltak", "epsilonk", "site", "Deltai", "check");
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    epsilonkprime = - t1 * cos(k) - t2 * cos(2.0 * k) - mu;
    Deltai = 0.0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      Deltakprime   = gsl_vector_get(D, iprime);
      Deltai += 1.0 / (Ndouble) * (Deltakprime * sin(kprime * i));
    }
    printf("%lg \t %lg \t %lg \t %lg \t %lg \t %i \n", k, gsl_vector_get(D, i), epsilonkprime, (double)i, Deltai, check);
  }
  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  return 0;
}