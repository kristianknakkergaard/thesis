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
double VFF0(double k, double rBB, double rBF, double nB, double mB, double lt)
{
  double xi     = M_PI / 2.0 * sqrt(1.0 / rBB ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/ 

  double factor = 4.0 / M_PI * (mB + 1.0/mB + 2.0) * 1.0 / (lt * lt) * rBF * rBF * pow(nB, 1.0 / 3.0) ;
  double f      = - factor * 1.0 / (k * k + 2.0 / (xi * xi) );
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
  double lt  = 0.5;

  /*nB:*/
  double nB  = 10.0; 

  /*k-values:*/
  double k_low = -40.0, k_up = 40.0, dk = 0.01;
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
  gsl_matrix *VFF0matrix = gsl_matrix_calloc(N, N);

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
      gsl_matrix_set( VFF0matrix, i, iprime, VFF0(k - kprime, rBB, rBF, nB, mB, lt) );
    }
  }

  /*For calculating mu:*/
  double mu_low = 0.7, dmu = 0.0003; 
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

  /*k for maximum of Delta*/
  int kmax_index;
  double kmax;

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
        epsilonkprime = kprime * kprime - mu;
        Deltakprime   = gsl_vector_get(D, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
        gapintegrand  =  -1.0 / (2.0 * M_PI) * gsl_matrix_get(VFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
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
        muintegrand   = 1.0 / 2.0 * 1.0/2.0 * (1.0 -  epsilonkprime / EFkprime);
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
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k, gsl_vector_get(D,i), mu);
  }
  printf("\n\n");

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

  /*Temperatures:*/
  const int numberT = 4; 
  double T[numberT] = {0.075, 0.100, 0.120, 0.127};

  for (int iT = 0; iT < numberT; ++iT)
  {
    check = 0;
    for (int j = 0; j < 1000; ++j)
    {
      D_maxk = gsl_vector_max(D);
      for (int i = 0; i < N; ++i)
      {
        if (D_maxk == 0)
        {
          break;
        }
        k = ((double)i) * dk + k_low;
        gapintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = kprime * kprime - mu;
          Deltakprime   = gsl_vector_get(D, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  =  -1.0 / (2.0 * M_PI) * gsl_matrix_get(VFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime) * tanh(EFkprime / (2.0 * T[iT]));
          gapintegral  += gapintegrand*dk; 
        }
        gsl_vector_set(D, i, gapintegral);
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
          muintegrand   = 1.0 / 2.0 * ( epsilonkprime / EFkprime * ( 1.0/(exp(EFkprime / T[iT]) + 1.0) - 1.0/2.0) + 1.0/2.0 );
          muintegral   += muintegrand * dk;
        }
        gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
      {
        break;
      }
      mu = mu_min_value;
      ++check;
    }

    kmax_index  = gsl_vector_max_index(D);
    kmax        = ((double)kmax_index) * dk + k_low;
    fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %i \n", T[iT], kmax, gsl_vector_max(D), mu_min_value, check);

    
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg  \n", k, gsl_vector_get(D,i), T[iT], mu_min_value );
    }

    printf("\n\n");
  }

  
  gsl_vector_free(D);
  gsl_matrix_free(VFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}