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
  double f      = - factor * gsl_sf_bessel_K0(sqrt(pow(q * d, 2.0) + 2.0 * pow( d/xi, 2.0)  ));
  return f;
}

double Deltaguess11(double k)
{
  return 0.4 * k / (pow(k, 4.0) + 1);
}

double Deltaguess12(double k)
{
  return 0.3 * exp( - 0.2 * k * k);
}


int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.11;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/
  double d_low = 0.7, d_up = 1.0, dd = 0.01; 

  double convergence = 5e-4;

  /*nB:*/
  double nB  = 100.0;

  /*k-values:*/
  double k_low = -25.0, k_up = 25.0, dk = 0.05;
  int N = (int) (k_up - k_low)/dk;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double EFminuskprime;
  double Delta11kprime;
  double Delta12kprimereal;
  double Delta12kprimeimag;

  double gapintegrand11;
  double gapintegrand12real;
  double gapintegrand12imag;
  double gapintegrandcheck;

  double gapintegral11;
  double gapintegral12real;
  double gapintegral12imag;
  double gapintegralcheck;

  double D11_maxk;
  double D12real_maxk;
  double D12imag_maxk;
  double D12ratio;

  /*vectors and matrices:*/
  gsl_vector *D11     = gsl_vector_calloc(N);
  gsl_vector *D12real = gsl_vector_calloc(N);
  gsl_vector *D12imag = gsl_vector_calloc(N);
  gsl_vector *Dcheck  = gsl_vector_calloc(N);
  gsl_vector *E       = gsl_vector_calloc(N);

  gsl_matrix *WFF11matrix = gsl_matrix_calloc(N, N);
  gsl_matrix *VFF12matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;
  int check = 0;
  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF11matrix, i, iprime, WFF11(k,kprime, rBB, rBF, nB, mB) );
    }
  }

  /*For calculating mu:*/
  double mu_low = 0.9, dmu = 0.0001;
  double muintegral, muintegrand;
  double mu_min_value;
  double mu = 1.0;
  int Nmu = 3000;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*T = 0: */
  printf("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n", "k", "Delta11", "Delta12real", "Delta12imag", "Emin", "mu", "iterations", "d", "D12ratio");

  printf("\n\n");
  for (double d = d_low; d < d_up; d+=dd)
  {
    check = 0;
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set(D11, i, Deltaguess11(k));
      gsl_vector_set(D12real, i, 0.2 * Deltaguess12(k));
      gsl_vector_set(D12imag, i, Deltaguess12(k));

      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( VFF12matrix, i, iprime, VFF12(k - kprime, d, rBB, rBF, nB, mB) );
      }
    }
    for (int j = 0; j < 10000; ++j)
    {
      D11_maxk     = gsl_vector_max(D11);
      D12real_maxk = gsl_vector_max(D12real);
      D12imag_maxk = gsl_vector_max(D12imag);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral11     = 0.0;
        gapintegral12real = 0.0;
        gapintegral12imag = 0.0;
        gapintegralcheck  = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime              = ((double)iprime) * dk + k_low;
          epsilonkprime       = kprime * kprime - mu;
          Delta11kprime       = gsl_vector_get(D11, iprime);
          Delta12kprimereal   = gsl_vector_get(D12real, iprime);
          Delta12kprimeimag   = gsl_vector_get(D12imag, iprime);
          EFkprime            = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + Delta12kprimeimag * Delta12kprimeimag  - 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
          EFminuskprime            = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + Delta12kprimeimag * Delta12kprimeimag  + 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
          
          gapintegrand11      =  -1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF11matrix, i, iprime) * (Delta11kprime - Delta12kprimereal) / (2.0 * EFkprime);
          gapintegrand12real  =  -1.0 / (2.0 * M_PI) * gsl_matrix_get(VFF12matrix, i, iprime) * ( - (Delta11kprime - Delta12kprimereal) / (4.0 * EFkprime) + (Delta11kprime + Delta12kprimereal) / (4.0 * EFminuskprime) );

          gapintegrand12imag  = -1.0 / (2.0 * M_PI) * gsl_matrix_get(VFF12matrix, i, iprime) * Delta12kprimeimag * ( 1.0/ (4.0 * EFkprime) + 1.0 / (4.0 * EFminuskprime) );

          gapintegrandcheck   = 1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF11matrix, i, iprime) * Delta12kprimeimag / (2.0 * EFkprime);

          gapintegral11      += gapintegrand11*dk; 
          gapintegral12real  += gapintegrand12real*dk;
          gapintegral12imag  += gapintegrand12imag*dk;
          gapintegralcheck   += gapintegrandcheck*dk; 
        }
        gsl_vector_set(D11,    i, gapintegral11);
        gsl_vector_set(D12real,i, gapintegral12real);
        gsl_vector_set(D12imag,i, gapintegral12imag);
        gsl_vector_set(Dcheck, i, gapintegralcheck);
      }
      for (int imu = 0; imu < Nmu; ++imu)
      {
        mu_guess = (double) mu_low + imu * dmu;
        muintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime              = ((double)iprime) * dk + k_low;
          epsilonkprime       = kprime * kprime - mu_guess;
          Delta11kprime       = gsl_vector_get(D11, iprime);
          Delta12kprimereal   = gsl_vector_get(D12real, iprime);
          Delta12kprimeimag   = gsl_vector_get(D12imag, iprime);
          EFkprime            = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + Delta12kprimeimag * Delta12kprimeimag  - 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);

          muintegrand     = 1.0/2.0 * 1.0 / 2.0 * (1.0 -  epsilonkprime / EFkprime);
          muintegral     += muintegrand*dk;
        }
        gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
      }
      mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence && fabs(gsl_vector_max(D12real) - D12real_maxk)/D12real_maxk < convergence && fabs(gsl_vector_max(D12imag) - D12imag_maxk)/D12imag_maxk < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence && gsl_vector_max(D12real) < convergence && fabs(gsl_vector_max(D12imag) - D12imag_maxk)/D12imag_maxk < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence && gsl_vector_max(D12imag) < convergence && fabs(gsl_vector_max(D12real) - D12real_maxk)/D12real_maxk < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence && gsl_vector_max(D12imag) < convergence && gsl_vector_max(D12real) < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence && fabs(gsl_vector_max(D12real) - D12real_maxk)/D12real_maxk < convergence && fabs(gsl_vector_max(D12imag) - D12imag_maxk)/D12imag_maxk < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence && gsl_vector_max(D12real) < convergence && fabs(gsl_vector_max(D12imag) - D12imag_maxk)/D12imag_maxk < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence && gsl_vector_max(D12imag) < convergence && fabs(gsl_vector_max(D12real) - D12real_maxk)/D12real_maxk < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence && gsl_vector_max(D12real) < convergence && gsl_vector_max(D12imag) < convergence && gsl_vector_max(Dcheck) < convergence)
      {
        break;
      }
      mu = mu_min_value;
      ++check; 
    }
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime              = ((double)iprime) * dk + k_low;
      epsilonkprime       = kprime * kprime - mu;
      Delta11kprime       = gsl_vector_get(D11, iprime);
      Delta12kprimereal   = gsl_vector_get(D12real, iprime);
      Delta12kprimeimag   = gsl_vector_get(D12imag, iprime);
      EFkprime            = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + Delta12kprimeimag * Delta12kprimeimag  - 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
      
      gsl_vector_set(E, iprime, EFkprime);
    }

    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      if (gsl_vector_max(D12imag) > convergence )
      {
        D12ratio = gsl_vector_max(D12real)/gsl_vector_max(D12imag);
        printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %i \t %lg \t %lg \n", k, gsl_vector_get(D11,i), gsl_vector_get(D12real,i), gsl_vector_get(D12imag,i), gsl_vector_min(E), mu, check, d, D12ratio);
      }
      else
        printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %i \t %lg \t %lg \n", k, gsl_vector_get(D11,i), gsl_vector_get(D12real,i), gsl_vector_get(D12imag,i), gsl_vector_min(E), mu, check, d, 0.0);
    }
    printf("\n \n");
  }
  


  
  gsl_vector_free(D11);
  gsl_vector_free(D12real);
  gsl_vector_free(D12imag);
  gsl_vector_free(Dcheck);
  gsl_vector_free(E);
  gsl_matrix_free(WFF11matrix);
  gsl_matrix_free(VFF12matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}