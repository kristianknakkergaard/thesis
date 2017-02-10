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

double WFF12(double k, double kprime, double d, double rBB, double rBF, double nB, double mB)
{

  return 1.0 / 2.0 * ( VFF12(k + kprime, d, rBB, rBF, nB, mB) + VFF12(k - kprime, d, rBB, rBF, nB, mB) ); 
}

double Deltaguess11(double k)
{
  return 0.2 * k / (pow(k, 4.0) + 1);
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
  double d_low = 0.74, d_up = 0.75, dd = 0.0002;  

  double convergence = 1e-4;

  double max_iter = 1000;

  double nB  = 100.0; 

  /*k-values:*/
  double k_low = -25.0, k_up = 25.0, dk = 0.05;
  int N = (int) (k_up - k_low)/dk;

  /*variables:*/
  double epsilonkprime;
  double EFkprimeminus;
  double EFkprimeplus;
  double Delta11kprime;
  double Delta12kprimereal;

  double gapintegrand11;
  double gapintegrand12real;

  double gapintegral11;
  double gapintegral12real;

  double D11_maxk;
  double D12real_maxk;

  double E0integrand;
  double E0integral2;

  /*vectors and matrices:*/
  gsl_vector *D11     = gsl_vector_calloc(N);
  gsl_vector *D12real = gsl_vector_calloc(N);
  
  gsl_matrix *WFF11matrix = gsl_matrix_calloc(N, N);
  gsl_matrix *WFF12matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;
  int check2 = 0;
  int d_iter = 0;

  printf("%s \t %s \t %s \t %s \t %s \t %s \n", "d", "E02", "mu2", "Delta11max", "Delta12max", "check2");

  printf("\n\n");

  /*For calculating mu:*/
  double mu_low = 0.95, dmu = 0.00006;
  double muintegral, muintegrand;
  double mu_min_value;
  double mu = 1.0;
  double mu2; 
  int Nmu = 3000;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /* d-independent interaction: */
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF11matrix, i, iprime, WFF11(k, kprime, rBB, rBF, nB, mB) );
    }
  }

  for (double d = d_low; d < d_up; d+=dd)
  {
    check2 = 0;

    /*Initial guess for pairings*/
    /*Real interwire pairing*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set(D11,     i, 2.82 * Deltaguess11(k));
      gsl_vector_set(D12real, i, Deltaguess12(k));
      
      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF12matrix, i, iprime, WFF12(k, kprime, d, rBB, rBF, nB, mB) );
      }
    }

    mu = 1.0;

    for (int j = 0; j < max_iter; ++j)
    {
      D11_maxk     = gsl_vector_max(D11);
      D12real_maxk = gsl_vector_max(D12real);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral11     = 0.0;
        gapintegral12real = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime              = ((double)iprime) * dk + k_low;
          epsilonkprime       = kprime * kprime - mu;
          Delta11kprime       = gsl_vector_get(D11, iprime);
          Delta12kprimereal   = gsl_vector_get(D12real, iprime);

          EFkprimeminus       = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal - 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
          EFkprimeplus        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
          
          gapintegrand11      = -1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF11matrix, i, iprime) * ( (Delta11kprime + Delta12kprimereal) / (4.0 * EFkprimeplus) + (Delta11kprime - Delta12kprimereal) / (4.0 * EFkprimeminus) );
          gapintegrand12real  = -1.0 / (2.0 * M_PI) * gsl_matrix_get(WFF12matrix, i, iprime) * ( (Delta12kprimereal + Delta11kprime) / (4.0 * EFkprimeplus) + (Delta12kprimereal - Delta11kprime) / (4.0 * EFkprimeminus) );

          gapintegral11      += gapintegrand11*dk; 
          gapintegral12real  += gapintegrand12real*dk;
        }
        gsl_vector_set(D11,    i, gapintegral11);
        gsl_vector_set(D12real,i, gapintegral12real);
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

          EFkprimeminus       = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal - 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
          EFkprimeplus        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);

          muintegrand         =  1.0 / 8.0 * (2.0 -  epsilonkprime * ( 1.0 / EFkprimeminus + 1.0 / EFkprimeplus ) );
          muintegral         += muintegrand * dk;
        }
        gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence && fabs(gsl_vector_max(D12real) - D12real_maxk)/D12real_maxk < convergence )
      {
        break;
      }
      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D11) - D11_maxk)/D11_maxk < convergence && gsl_vector_max(D12real) < convergence )
      {
        break;
      }
      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence && fabs(gsl_vector_max(D12real) - D12real_maxk)/D12real_maxk < convergence )
      {
        break;
      }
      if ( fabs(mu - mu_min_value)/mu_min_value < convergence && gsl_vector_max(D11) < convergence && gsl_vector_max(D12real) < convergence )
      {
        break;
      }
      mu = mu_min_value;
      ++check2; 
    }
    mu2 = mu; 

    E0integral2 = 0.0;

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime              = ((double)iprime) * dk + k_low;
      epsilonkprime       = kprime * kprime - mu;
      Delta11kprime       = gsl_vector_get(D11, iprime);
      Delta12kprimereal   = gsl_vector_get(D12real, iprime);

      EFkprimeminus       = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal - 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);
      EFkprimeplus        = pow(epsilonkprime * epsilonkprime + Delta11kprime * Delta11kprime + Delta12kprimereal * Delta12kprimereal + 2.0 * Delta11kprime * Delta12kprimereal, 1.0 / 2.0);

      E0integrand         = - 1.0 / 8.0 * ( pow(epsilonkprime - EFkprimeplus, 2.0) / EFkprimeplus + pow(epsilonkprime - EFkprimeminus, 2.0) / EFkprimeminus ); 
      E0integral2        += dk * E0integrand;
    }

    printf("%lg \t %lg \t %lg \t %lg \t %lg \t %i \n", d, E0integral2 + 2.0 * mu2, mu2, gsl_vector_max(D11), gsl_vector_max(D12real), check2);
    fprintf(stderr, "d_iter = %i\n", ++d_iter );
  }


  gsl_vector_free(D11);
  gsl_vector_free(D12real);

  gsl_matrix_free(WFF11matrix);
  gsl_matrix_free(WFF12matrix);

  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}