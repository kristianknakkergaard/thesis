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
  double xi     = M_PI / sqrt( 8.0 * nB * aB ); /*xi * kF*/ 

  double factor = 4.0 / (M_PI * M_PI) * (1.0/mB + mB + 2.0) * nB * pow(aBF, 2.0);
  double f      = - factor * log( ( pow(k + kprime, 2.0) + 2.0/(xi * xi) )/( pow(k - kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}

double Deltaasymp(double D_maxkT, double T, double TC)
{
  return D_maxkT * pow(1.0 - pow(T / TC, 3.0) , 1.0/2.0);
}

double Deltaguess(double k)
{
  return 0.4 * k / ( pow(k, 4.0) + 1.0);
}


int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (3 percent depletion) */
  double rBF = 0.07;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  /*double inverseBdist = 0.036;
  double nB           = pow(inverseBdist, -3.0); */
  double nB  = 1.0e5;

  /*k-values:*/
  double k_low = 0.0, k_up = 400.0, dk = 0.8;
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
    gsl_vector_set(D, i, Deltaguess(k));

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, rBB, rBF, nB, mB) );
    }
  }

  /*For calculating mu:*/
  double mu_low = 0.4, dmu = 0.001;
  double muintegral, muintegrand;
  double mu_min_value;
  double mu = 1.0;
  int Nmu = 900;
  double mu_guess;
  gsl_vector*   mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector* find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*k for maximum of Delta*/
  int kmax_index;
  double kmax;

  /*T = 0: */
  for (int j = 0; j < 300; ++j)
  {
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

        gapintegral  += gapintegrand * dk; 
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

        muintegral   += muintegrand * dk;
      }
      gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
    }

    mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));

    if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
    {
      break;
    }
    mu     = mu_min_value;
    D_maxk = gsl_vector_max(D);
    kmax   = (double)gsl_vector_max_index(D) * dk + k_low;

    ++check; 
  }

  double D_maxkT  = D_maxk; 
  kmax            = (double)gsl_vector_max_index(D) * dk + k_low;
  double aB       = M_PI * rBB / pow(nB, 1.0/3.0);
  double aBF      = M_PI * rBF / pow(nB, 1.0/3.0);
  double retneg   = pow(1.0/mB, 2.0) * nB * aB;

  fprintf(stderr, "(nBaB^3)^(1/3) = %lg, (nBaBF^3)^(1/3) = %lg, mB/mF = %lg, nB/nF^3 = %lg, kF*aB = %lg, kF*aBF = %lg, (mF/mB)^2*nB/nF^3*kF*aB = %lg \n", rBB, rBF, mB, nB, aB, aBF, retneg);

  fprintf(stderr, "\n \n");

  fprintf(stderr, "%s \t %s \t %s \t %s \t %s \n", "T", "kmax", "Deltamax", "mu", "check");
  fprintf(stderr, "\n \n");

  fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %i \n", 0.0, kmax, gsl_vector_max(D), mu, check);

  /*We list the function values:*/
  for (int i = 0; i < N-1; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), 0.0);
  }
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k, gsl_vector_get(D,i), 0.0);
  }
  printf("\n\n");

  /*Temperatures:*/
  double T_low = 0.0001, T_high = 1.0, dT = 0.0001;

  double TC = 0.0;
  for (double T = T_low; T < T_high; T+=dT)
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
          gapintegrand  =  -1.0/M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime) * tanh(EFkprime / (2.0 * T));

          gapintegral  += gapintegrand * dk; 
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
          muintegrand   = epsilonkprime / EFkprime * ( 1.0/(exp(EFkprime / T) + 1.0) - 1.0/2.0) + 1.0/2.0;

          muintegral   += muintegrand * dk;
        }
        gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
      }
      mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs( gsl_vector_max(D) - D_maxk )/D_maxk < 1e-3 )
      {
        break;
      }
      mu = mu_min_value;
      ++check;
    }

    if (gsl_vector_max(D)/D_maxkT < 1e-2)
    {
      TC = T;
      break;
    }
    kmax_index  = gsl_vector_max_index(D);
    kmax        = ((double)kmax_index) * dk + k_low;

    fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %i \n", T, kmax, gsl_vector_max(D), mu_min_value, check);

    for (int i = 0; i < N-1; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), T, mu_min_value);
    }
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg  \n", k, gsl_vector_get(D,i), T, mu_min_value );
    }

    printf("\n\n");
  }
  fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %i \n", TC, 0.0, 0.0, mu_min_value, 0);

  /****************************************/
  fprintf(stderr, "\n\n");
  for (double T = 0; T < TC; T+=dT)
  {
    fprintf(stderr, "%lg \t %lg \n", T, Deltaasymp(D_maxkT, T, TC) ); 
  }
  fprintf(stderr, "%lg \t %lg \n", TC, 0.0);

  
  gsl_vector_free (D);
  gsl_matrix_free (WFF0matrix);
  gsl_vector_free (mu_vector);
  gsl_vector_free (find_min_mu);
  return 0;
}