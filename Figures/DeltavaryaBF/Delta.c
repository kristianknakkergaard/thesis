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
double WFF0(double k, double kprime, double BFscatlength)
{
  double densityratio = 1.0;   /*nB/nF^3 */
  double BBscatlength = 0.001;
  double xi           = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/ 
  double massratio    = 0.8; /*mB/mF*/
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
  double k_low = 0.0, k_up = 50.0, dk = 0.1;
  int N = (int) (k_up - k_low)/dk;
  /*fermion mass*/

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double gapintegrand;
  double gapintegral;

  /*BFscatlengths:*/
  double BFscatlength_low = 0.150, BFscatlength_up = 0.2, dBFscatlength = 0.0001;

  gsl_vector *D = gsl_vector_calloc(N);
  double D_maxk;
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;

  /*For calculating mu:*/
  double mu_low = 0.8, dmu = 0.001;
  double muintegral, muintegrand;
  double mu_min_value = 0.98;
  int Nmu = 250;
  double mu, mu_guess;
  int mu_check; 

  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  /*mu_vector*/
  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }


  /*Delta guess:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k));
  }

  /*Number of iterations*/
  int iter_number = 0;
  
  for (double BFscatlength = BFscatlength_up; BFscatlength > BFscatlength_low; BFscatlength-=dBFscatlength)
  {
    mu_check = 0;
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, BFscatlength) );
      }
    }
    for (int j = 0; j < 500; ++j)
    {
      mu = mu_min_value;
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
        gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));
      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
      {
        break;
      }
      ++mu_check; 
    }
    fprintf(stderr, "%i \t %i\n", ++iter_number, mu_check);
    printf("%lg \t %lg \t %lg \t %lg \n", BFscatlength, gsl_vector_max(D), mu_min_value, 0.0);
  }

  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);

  return 0;
}