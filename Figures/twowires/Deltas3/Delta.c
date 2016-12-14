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

/*nB/nF^3, kF * aB HERE*/
/*SET mB/mF, nB/nF^3 and kF * aBF */
double WFF11(double k, double kprime, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 

  double factor = 4.0/(M_PI*M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);
  double f      = - factor * log( ( pow(k+kprime, 2.0) + 2.0/(xi * xi) )/( pow(k-kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}

double VFF12(double q, double d, double rBB, double rBF, double nB, double mB)
{
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);    /*kF * aB*/
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0); /*kF * aBF */
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 

  double factor = 16.0/(M_PI * M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);
  double f      = - factor * gsl_sf_bessel_K0(sqrt(pow(q * d, 2.0) + pow(2.0 * d/xi, 2.0)  ));
  return f;
}

double Deltaasymp(double D11_maxkT, double T, double TC)
{
  return D11_maxkT * pow(1.0 - pow(T/TC, 3.0) , 1.0/2.0);
}

double Deltaguess11(double k)
{
  return 0.4 * k / (pow(k, 4.0) + 1);
}

double Deltaguess12(double k)
{
  return 0.4 * k / (pow(k, 4.0) + 1);
}


int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.11;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/
  double d_low = 0.001, d_up = 0.68, dd = 0.001; 

  /*interparticle distances::*/
  double Bdist_min = 0.01, Bdist_max = 0.8, dBdist = 0.005; 
  double nB = 0.0;

  /*k-values:*/
  double k_low = -25.0, k_up = 25.0, dk = 0.1;
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

  /*vectors and matrices:*/
  gsl_vector *D11 = gsl_vector_calloc(N);
  gsl_vector *D12 = gsl_vector_calloc(N);
  double D11_maxk;
  double D12_maxk;
  gsl_matrix *WFF11matrix = gsl_matrix_calloc(N, N);
  gsl_matrix *VFF12matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;
  int check = 0;
  /*start guess for Delta:*/
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

  /*For calculating mu:*/
  double mu_low = - 0.1, dmu = 0.001;
  double muintegral, muintegrand;
  double mu_min_value;
  double mu = 1.0;
  int Nmu = 1500;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*T = 0: */
  int    d_iter  = 0;
  int    nB_iter = 0;
  double d_equal = 0.0;

  printf("%s \t %s \t %s \t %s\n", "nF/nB^(1/3)", "d_equal", "Delta11max", "Delta12max");

  for (double Bdist = Bdist_max; Bdist > Bdist_min; Bdist -= dBdist)
  {
    nB = pow(Bdist, -3.0);
    
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF11matrix, i, iprime, WFF11(k,kprime, rBB, rBF, nB, mB) );
      }
    }

    d_iter = 0;

    for (double d = d_up; d > d_low; d -= dd)
    {

      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gsl_vector_set(D11, i, Deltaguess11(k));
        gsl_vector_set(D12, i, Deltaguess12(k));
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime = ((double)iprime) * dk + k_low;
          gsl_matrix_set( VFF12matrix, i, iprime, VFF12(k - kprime, d, rBB, rBF, nB, mB) );
        }
      }
      for (int j = 0; j < 1000; ++j)
      {
        D11_maxk = gsl_vector_max(D11);
        D12_maxk = gsl_vector_max(D12);
        for (int i = 0; i < N; ++i)
        {
          k = ((double)i) * dk + k_low;
          gapintegral11 = 0.0;
          gapintegral12 = 0.0;
          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime          = ((double)iprime) * dk + k_low;
            epsilonkprime   = kprime * kprime - mu;
            Delta11kprime   = gsl_vector_get(D11, iprime);
            Delta12kprime   = gsl_vector_get(D12, iprime);
            EFkprime        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprime * Delta12kprime, 1.0 / 2.0);
            gapintegrand11  =  -1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF11matrix, i, iprime) * Delta11kprime / (2.0 * EFkprime);
            gapintegrand12  =  -1.0 / (2.0 * M_PI) * gsl_matrix_get(VFF12matrix, i, iprime) * Delta12kprime / (2.0 * EFkprime);
            gapintegral11  += gapintegrand11*dk; 
            gapintegral12  += gapintegrand12*dk; 
          }
          gsl_vector_set(D11,i,gapintegral11);
          gsl_vector_set(D12,i,gapintegral12);
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
            Delta12kprime   = gsl_vector_get(D12, iprime);
            EFkprime        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprime * Delta12kprime, 1.0 / 2.0);
            muintegrand     = 1.0/2.0 * 1.0 / 2.0 * (1.0 -  epsilonkprime / EFkprime);
            muintegral     += muintegrand*dk;
          }
          gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
        }
        mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

        if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < 1e-3 && fabs(gsl_vector_max(D12) - D12_maxk)/D12_maxk < 1e-3)
        {
          break;
        }

        if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < 1e-3 && gsl_vector_max(D12) < 1e-3)
        {
          break;
        }

        if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && gsl_vector_max(D11) < 1e-3 && fabs(gsl_vector_max(D12) - D12_maxk)/D12_maxk < 1e-3)
        {
          break;
        }

        if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && gsl_vector_max(D11) < 1e-3 && gsl_vector_max(D12) < 1e-3)
        {
          break;
        }
        mu = mu_min_value;
        ++check; 
      }
      if (gsl_vector_max(D11) < gsl_vector_max(D12))
      {
        d_equal = d;
        d_up    = d + dd;
        break;
      }
      fprintf(stderr, "d_iter = %i \t nB_iter = %i \n", ++d_iter, nB_iter );
    }
    printf("%lg \t %lg \t %lg \t %lg \n", Bdist, d_equal, gsl_vector_max(D11), gsl_vector_max(D12));
    ++nB_iter;
  }
  
  
  gsl_vector_free(D11);
  gsl_vector_free(D12);
  gsl_matrix_free(WFF11matrix);
  gsl_matrix_free(VFF12matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}