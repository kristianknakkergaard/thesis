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
double WFF0(double k, double kprime)
{
  double BBscatlength = 0.3; /*kF * aB*/
  double densityratio = 1.0;   /*nB/nF^3 */
  double xi           = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/ 
  double massratio    = 0.8; /*mB/mF*/
  double BFscatlength = 0.3; /*kF * aBF */
  double factor       = pow(2.0,5.0)/(M_PI*M_PI) * pow(BFscatlength, 2.0) * densityratio * (massratio + 1.0/massratio + 2.0);

  double f = - factor * log( ( pow(k+kprime, 2.0) + 2.0/(xi * xi) )/( pow(k-kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}


double Deltaguess(double k)
{
  return 0.4 * k / (pow(k,4) + 1);
}


int
main (void)
{
  /*k-values:*/
  double k_low = 0.0, k_up = 70.0, dk = 0.01;
  int N = (int) (k_up - k_low)/dk;
  /*fermion mass*/

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double gapintegrand;
  double gapintegral;
  double mu = 1.0;

  gsl_vector *D = gsl_vector_calloc(N);
  double D_max;
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;


  /*For calculating mu:*/
  double mu_low = 0.5, mu_up = 1.5, dmu = 0.001;
  double muintegral, muintegrand;
  double mu_min_value;
  int Nmu = (int) (mu_up - mu_low)/dmu;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k));

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k,kprime) );
    }
  }
  
  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }
  
  /*T = 0: */
  for (int iter = 0; iter < 30; ++iter)
  {
    D_max = gsl_vector_max(D);
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
      gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
    }
    mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

    if ( fabs(mu - mu_min_value)/mu_min_value < 1e-2 && fabs(gsl_vector_max(D) - D_max)/D_max < 1e-2 )
    {
      break;
    }
    mu = mu_min_value;
  }

  double kmax = (double)gsl_vector_max_index(D) * dk + k_low;
  fprintf(stderr, "%lg \t %lg \t %lg \t %lg \n", 0.0, kmax, gsl_vector_max(D), mu);


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

  /*Now: temperature T>0 */
  double T_low = 0.002, T_high = 0.2, dT = 0.002;
  int kmax_index;

  for (double T = T_low; T < T_high; T+=dT)
  {
    for (int iter = 0; iter < 30; ++iter)
    {
      D_max = gsl_vector_max(D);
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
          gapintegrand  =  -1.0/M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0*EFkprime) * tanh(EFkprime / (2.0*T));
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
          muintegrand   = epsilonkprime / EFkprime * ( 1.0/(exp(EFkprime / T) + 1.0) - 1.0/2.0) + 1.0/2.0;
          muintegral   += muintegrand*dk;
        }
        gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-2 && fabs(gsl_vector_max(D) - D_max)/D_max < 1e-2 )
      {
        break;
      }
      mu = mu_min_value;
    }

    kmax_index  = gsl_vector_max_index(D);
    kmax        = ((double)kmax_index) * dk + k_low;

    fprintf(stderr, "%lg \t %lg \t %lg \t %lg \n", T, kmax, gsl_vector_max(D), mu);
    for (int i = 0; i < N-1; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), T, mu);
    }
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg  \n", k, gsl_vector_get(D,i), T, mu );
    }

    printf("\n\n");
  }

  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}