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
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/ 

  double factor = 4.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;
  double f      = - factor * log( ( pow( k + kprime, 2.0 ) + 2.0 / (xi * xi) ) / ( pow( k - kprime, 2.0 ) + 2.0 / (xi * xi) ) );
  return f; 
}

double Deltaasymp(double D_maxkT, double T, double TC)
{
  return D_maxkT * pow(1.0 - pow(T/TC, 3.0) , 1.0/2.0);
}

double Deltaguess(double k, double rBB, double rBF, double nB, double mB)
{
  return k / (k * k + 3.3147);
}


int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  double nB  = 100.0; 

  /*k-values:*/
  double k_low = 0.0, k_up = 25.0, dk = 0.01;
  int N = (int) (k_up - k_low)/dk;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double gapintegrand;
  double gapintegral;

  /*vectors and matrices:*/
  gsl_vector *D = gsl_vector_calloc(N);
  double D_maxk;
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;
  int check = 0;
  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k, rBB, rBF, nB, mB));

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k,kprime, rBB, rBF, nB, mB) );
    }
  }

  /*For calculating mu:*/
  double mu_low = 0.95, dmu = 0.00005;
  double muintegral, muintegrand;
  double mu_min_value;
  double mu = 1.0;
  int Nmu = 2000;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*T = 0: */
  for (int j = 0; j < 50; ++j)
  {
    D_maxk = gsl_vector_max(D);
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gapintegral = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = kprime * kprime - mu;
        Deltakprime   = gsl_vector_get(D, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
        gapintegrand  =  -1.0 / M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
        gapintegral  += gapintegrand*dk; 
      }
      gsl_vector_set(D,i,gapintegral);
    }
    for (int imu = 0; imu < Nmu; ++imu)
    {
      mu_guess = (double) mu_low + imu * dmu;
      muintegral = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = kprime * kprime - mu_guess;
        Deltakprime   = gsl_vector_get(D, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
        muintegrand   = 1.0/2.0 * (1.0 -  epsilonkprime / EFkprime);
        muintegral   += muintegrand*dk;
      }
      gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
    }
    mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

    if ( fabs(mu - mu_min_value)/mu_min_value < 1e-4 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-4 )
    {
      break;
    }
    mu = mu_min_value;
    ++check; 
  }

  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0);
  double retneg = pow(1.0/mB, 2.0) * nB * aB;
  fprintf(stderr, "(nBaB^3)^(1/3) = %lg, (nBaBF^3)^(1/3) = %lg, mB/mF = %lg, nB/nF^3 = %lg, kF*aB = %lg, kF*aBF = %lg, (mF/mB)^2*nB/nF^3*kF*aB = %lg \n", rBB, rBF, mB, nB, aB, aBF, retneg);

  fprintf(stderr, "\n \n");

  /*We list the function values:*/
  printf("%s \t %s \t %s \n", "k", "Deltak", "mu" );
  for (int i = 0; i < N-1; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), mu);
  }
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k, gsl_vector_get(D,i), mu);
  }
  /*Ground state energy: */

  double E0integrand; 
  double E0integral = 0; 
  for (int iprime = 0; iprime < N; ++iprime)
  {
    kprime        = ((double)iprime) * dk + k_low;
    epsilonkprime = kprime * kprime - mu;
    Deltakprime   = gsl_vector_get(D, iprime);
    EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

    E0integrand   = - 1.0 / 4.0 * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
    E0integral   += E0integrand*dk; 
  }

  fprintf(stderr, "E0 = %lg, mu = %lg, check = %i\n", E0integral + mu, mu, check);



  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}