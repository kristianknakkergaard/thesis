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
  double f      = - factor * gsl_sf_bessel_K0( sqrt( pow( q * d, 2.0 ) + 2.0 * pow( d / xi, 2.0 )  ) );
  return f;
}

double WFF12(double k, double kprime, double d, double rBB, double rBF, double nB, double mB)
{
  return 1.0 / 2.0 * ( VFF12(k + kprime, d, rBB, rBF, nB, mB) + VFF12(k - kprime, d, rBB, rBF, nB, mB) ); 
}

double Deltaguess11(double k)
{
  return 0.4 * k / (pow(k, 4.0) + 1);
}

double Deltaguess12(double k)
{
  return 0.4 / (pow(k, 4.0) + 1);
}


int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.11;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/
  double d_low = 0.001, d_up = 0.792, dd = 0.0005; 

  double convergence = 1e-4; 

  /*interparticle distances::*/
  double Bdist_min = 0.001, Bdist_max = 0.65, dBdist = 0.001; 
  double nB = 0.0;

  /*k-values:*/
  double k_low = -25.0, k_up = 25.0, dk = 0.05;
  int N = (int) (k_up - k_low)/dk;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Delta11kprime;
  double Delta12kprime;
  double gapintegrand11;
  double gapintegrand12;
  double gapintegral11;
  double gapintegral12;

  /*Minimal energy*/
  double E0integrand;
  double E0integral;
  double E01, E02; 

  /*vectors and matrices:*/
  gsl_vector *D11 = gsl_vector_calloc(N);
  gsl_vector *D12 = gsl_vector_calloc(N);

  double D11_maxk;
  double D12_maxk;
  
  gsl_matrix *WFF11matrix = gsl_matrix_calloc(N, N);
  gsl_matrix *WFF12matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;
  int check1 = 0, check2 = 0;

  /*For calculating mu:*/
  double mu_low = 0.95, dmu = 0.00006;
  double mu1 = 1.0, mu2 = 1.0;
  
  double mu_guess, muintegral, muintegrand, mu_min_value;
  
  int Nmu = 3000;

  gsl_vector *mu_vector   = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*T = 0: */
  int    d_iter  = 0;
  int    nB_iter = 0;
  double d_c     = 0.0;

  printf("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n", "nF/nB^(1/3)", "d_c", "E01", "E02", "Delta11max", "Delta12max", "check1", "check2");

  for (double Bdist = Bdist_max; Bdist > Bdist_min; Bdist -= dBdist)
  {
    nB = pow(Bdist, -3.0);
    
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set(D11, i, Deltaguess11(k));
      gsl_vector_set(D12, i, Deltaguess12(k));
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF11matrix, i, iprime, WFF11(k,kprime, rBB, rBF, nB, mB) );
      }
    }

    for (int j = 0; j < 1000; ++j)
    {
      D11_maxk = gsl_vector_max(D11);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral11 = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime          = ((double)iprime) * dk + k_low;
          
          epsilonkprime   = kprime * kprime - mu1;
          Delta11kprime   = gsl_vector_get(D11, iprime);
          EFkprime        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime, 1.0 / 2.0);

          gapintegrand11  = -1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF11matrix, i, iprime) * Delta11kprime / (2.0 * EFkprime);
          gapintegral11  += gapintegrand11*dk; 
        }
        gsl_vector_set(D11,i,gapintegral11);
      }
      for (int imu = 0; imu < Nmu; ++imu)
      {
        mu_guess = (double) mu_low + imu * dmu;
        muintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime          = ((double)iprime) * dk + k_low;
          
          epsilonkprime   = kprime * kprime - mu_guess;
          Delta11kprime   = gsl_vector_get(D11, iprime);
          EFkprime        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime, 1.0 / 2.0);
          
          muintegrand     = 1.0 / 4.0 * (1.0 -  epsilonkprime / EFkprime);
          muintegral     += muintegrand*dk;
        }
        gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu1 - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence)
      {
        break;
      }
      if ( fabs(mu1 - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence)
      {
        break;
      }
      mu1 = mu_min_value;
      ++check1; 
    }

    E0integral = 0.0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = kprime * kprime - mu1;
      Delta11kprime = gsl_vector_get(D11, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime, 1.0 / 2.0);
      
      E0integrand = -1.0 / 4.0 * pow(epsilonkprime - EFkprime, 2.0) / EFkprime; 
      E0integral += dk * E0integrand;
    }
    E01 = E0integral + 2.0 * mu1; 

    d_iter = 0;
    for (double d = d_up; d > d_low; d -= dd)
    {
      check2 = 0;
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime = ((double)iprime) * dk + k_low;
          gsl_matrix_set( WFF12matrix, i, iprime, WFF12(k, kprime, d, rBB, rBF, nB, mB) );
        }
      }
      for (int j = 0; j < 1000; ++j)
      {
        D12_maxk = gsl_vector_max(D12);
        for (int i = 0; i < N; ++i)
        {
          k = ((double)i) * dk + k_low;
          gapintegral12 = 0.0;
          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime          = ((double)iprime) * dk + k_low;
            epsilonkprime   = kprime * kprime - mu2;
            Delta12kprime   = gsl_vector_get(D12, iprime);
            EFkprime        = pow(epsilonkprime * epsilonkprime + Delta12kprime * Delta12kprime, 1.0 / 2.0);

            gapintegrand12  = -1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF12matrix, i, iprime) * Delta12kprime / (2.0 * EFkprime);
            gapintegral12  += gapintegrand12*dk; 
          }
          gsl_vector_set(D12, i, gapintegral12);
        }
        for (int imu = 0; imu < Nmu; ++imu)
        {
          mu_guess = (double) mu_low + imu * dmu;
          muintegral = 0.0;

          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime          = ((double)iprime) * dk + k_low;
            epsilonkprime   = kprime * kprime - mu_guess;
            Delta12kprime   = gsl_vector_get(D12, iprime);
            EFkprime        = pow(epsilonkprime * epsilonkprime + Delta12kprime * Delta12kprime, 1.0 / 2.0);
            muintegrand     = 1.0 / 4.0 * (1.0 -  epsilonkprime / EFkprime);
            muintegral     += muintegrand*dk;
          }
          gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
        }
        mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

        if ( fabs(mu2 - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D12) - D12_maxk)/D12_maxk < convergence)
        {
          break;
        }
        if ( fabs(mu2 - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D12) < convergence)
        {
          break;
        }
        mu2 = mu_min_value;
        ++check2; 
      }

      E0integral = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = kprime * kprime - mu2;
        Delta12kprime = gsl_vector_get(D12, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Delta12kprime * Delta12kprime, 1.0 / 2.0);
        
        E0integrand = -1.0 / 4.0 * pow(epsilonkprime - EFkprime, 2.0) / EFkprime; 
        E0integral += dk * E0integrand;
      }
      E02 = E0integral + 2.0 * mu2; 

      if (E02 < E01)
      {
        d_c  = d;
        d_up = d + 2.0 * dd;
        break;
      }
      fprintf(stderr, "d_iter = %i \t nB_iter = %i \n", ++d_iter, nB_iter );
    }
    printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %i \t %i \n", Bdist, d_c, E01, E02, D11_maxk, D12_maxk, check1, check2);
    ++nB_iter;
  }
  
  
  gsl_vector_free(D11);
  gsl_vector_free(D12);
  
  gsl_matrix_free(WFF11matrix);
  gsl_matrix_free(WFF12matrix);

  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}