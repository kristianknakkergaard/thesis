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
  double xi     = M_PI/sqrt(8.0 * nB * aB ); /*xi * kF*/ 

  double factor = 4.0/(M_PI*M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);
  double f      = - factor * log( ( pow(k+kprime, 2.0) + 2.0/(xi * xi) )/( pow(k-kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}


double Deltaguess(double k)
{
  return 0.4 * k / (pow(k,4) + 1);
}


int
main (void)
{
  /*variables */
  
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  double nB  = 1.0e5; 
  double aB, retneg;

  
  /*k-values:*/
  double k_low = 0.0, k_up = 150.0, dk = 0.005;
  int N = (int) (k_up - k_low)/dk;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double gapintegrand;
  double gapintegral;

  /*vectors and matrices:*/
  gsl_vector *D = gsl_vector_calloc(N);
  double kmax, D_maxk;
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);


  /*For calculating mu:*/
  double mu_low = 0.1, dmu = 0.001;
  double muintegral, vkprimenorm2;
  double mu_min_value = 0.9;
  double mu;
  int Nmu   = 1000;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  double k;
  double kprime;
  int check = 0;
  /*Occupancy plot */
  int N_plot = 20.0/dk;

  printf("%s \t %s \t %s \t %s \t %s \n", "k", "Occupancy", "mu", "kmax", "Energy" );

  for (int i = 0; i < N_plot; ++i)
  {
    k = ((double)i) * dk + k_low;
    if (k <= 1.0 )
    {
      printf("%lg \t %lg \t %lg \t %lg \t %lg \n", k, 1.0, 1.0, 0.0, fabs(k * k - 1.0));
    }
    else
    printf("%lg \t %lg \t %lg \t %lg \t %lg \n", k, 0.0, 1.0, 0.0, fabs(k * k - 1.0));
  }
  printf("\n \n");

  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k));
  }


  fprintf(stderr, "%i\n", ++check);

  /*gas parameter*/
  const int Nr = 4;
  int iter;
  double rBB[Nr] = {0.005, 0.006, 0.007, 0.01}; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */

  for (int r = 0; r < Nr; ++r)
  {
    aB     = M_PI * rBB[r] / pow(nB, 1.0/3.0);
    retneg = pow(1.0/mB, 2.0) * nB * aB;
    
    if (retneg < 20.0)
    {
      fprintf(stderr, "retneg = %lg \n", retneg );
      break;
    }
    
    /*start guess for Delta:*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set(D, i, Deltaguess(k));

      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, rBB[r], rBF, nB, mB) );
      }
    }

    iter = 0;
    for (int j = 0; j < 100; ++j)
    {
      D_maxk = gsl_vector_max(D);
      kmax   =  (double)gsl_vector_max_index(D) * dk;
      mu     = mu_min_value;

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
          gapintegrand  =  - 1.0 / M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          
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
          vkprimenorm2  = 1.0/2.0 * (1.0 -  epsilonkprime / EFkprime);

          muintegral   += vkprimenorm2 * dk;
        }
        gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs( gsl_vector_max(D) - D_maxk )/D_maxk < 1e-3 )
      {
        break;
      }
      fprintf(stderr, "%i \t %i\n", check, ++iter);
    }

    if (gsl_vector_max(D) < 1e-2)
    {
      fprintf(stderr, "%s\n","Pairing too small" );
      return 0;
    }
    
    for (int iprime = 0; iprime < N_plot; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = kprime * kprime - mu;
      Deltakprime   = gsl_vector_get(D, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

      vkprimenorm2  = 1.0 / 2.0 * ( 1.0 - epsilonkprime / EFkprime ); 
      printf("%lg \t %lg \t %lg \t %lg \t %lg \n", kprime, vkprimenorm2, mu, kmax, EFkprime);
    }
    printf("\n \n");
    fprintf(stderr, "%i \t %i\n", ++check, iter);
  }



  
  gsl_vector_free (D);
  gsl_matrix_free (WFF0matrix);
  gsl_vector_free (mu_vector);
  gsl_vector_free (find_min_mu);
  return 0;
}