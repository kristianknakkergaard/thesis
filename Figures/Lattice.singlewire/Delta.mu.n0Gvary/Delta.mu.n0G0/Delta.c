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


int
main (void)
{
  /*restricted variables */
  double t1       = 1.0;  /*the unit of energy */
  
  /*free parameters:*/
  double fill_low = 0.0, dfill = 0.005;
  int    N        = 300; 
  double Ndouble  = (double)N;
  int    Nfill    = 200;  
  
  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;
  
  /*function variables:*/
  double epsilonkprime;
  double kprime;
  int check = 0;

  /*for the chemical potential*/
  double mu_low = -3.0, mu, dmu = 0.001;
  int Nmu       = 6000;
  double muintegral, fillingnum = 0.0;
  double mu_min_value;
  
  double mu_guess;
  gsl_vector *mu_vector     = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu   = gsl_vector_calloc(Nmu);
  gsl_vector *fillingvector = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double)imu * dmu + mu_low;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }
 
  /* start of fill loop*/
  double fill;
  int checkfill = 0; 

  printf("%s \t %s \t %s \n", "fill", "mu", "muteo");

  for (int jfill = 0; jfill < Nfill; ++jfill)
  {
    fill = fill_low + (double)jfill * dfill; 
    check = 0;

    for (int imu = 0; imu < Nmu; ++imu)
    {
      mu_guess =  (double)imu * dmu + mu_low;
      muintegral = 0.0;

      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = - t1 * cos(kprime) - mu_guess;
        if (epsilonkprime < 0)
        {
          muintegral += 1.0 / Ndouble;
        }
      }
      gsl_vector_set(find_min_mu, imu, pow(fill - muintegral, 2.0) );
      gsl_vector_set(fillingvector, imu, muintegral );
    }
    mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));
    fillingnum = gsl_vector_get(fillingvector, gsl_vector_min_index(find_min_mu));

    mu = mu_min_value;
    ++check; 
    printf("%lg \t %lg \t %lg \n", fill, mu, -t1 * cos(M_PI * fill));
    fprintf(stderr, "checkfill = %i\n", ++checkfill );
  }
 
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  gsl_vector_free(fillingvector);
  return 0;
}