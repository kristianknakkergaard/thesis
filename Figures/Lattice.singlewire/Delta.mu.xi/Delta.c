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

double Vindx (double x, double rBB, double rBF, double nB, double mB, double eps0, double lt)
{
  double xi     = lt / 2.0 * sqrt(1.0 / rBB ); /*xi / a*/
  double factor = sqrt(2.0) * (mB + 1.0 / mB + 2.0) * rBF * rBF * 1.0 / (lt * lt * lt) * eps0 * 1.0 / nB * pow(rBB, -1.0 / 2.0);

  return - factor * exp(- sqrt(2.0) * fabs(x) / xi ); 
}

double WFF0(double k, double kprime, double rBB, double rBF, double nB, double mB, int N, double eps0, double lt)
{
  double W = 0.0;

  for (int i = 1; i <= N; ++i)
  {
    double x  = (double) i;
    W += 1.0 / 2.0 * ( cos((kprime - k) * i ) - cos( (kprime + k) * i ) ) * Vindx( x, rBB, rBF, nB, mB, eps0, lt); 
  }

  return W; 
}

double Deltaguess(double k)
{
  return sin(k);
}


int
main (void)
{
  /*restricted variables */
  double t1       = 1.0;  /*the unit of energy */
  double strength = 7.0; /*total strength of interaction*/
  double mB       = 7.0/40.0; /*mB/mF*/
  double rBF      = 0.005;    /*(nB * aBF^3)^(1/3) */
  double eps0     = 0.5; /*energy of particle box */
  double lt       = 0.1;
  double constant = sqrt(2.0) * (mB + 1.0 / mB + 2.0) * rBF * rBF * 1.0 / (lt * lt * lt) * eps0;  

  /*free parameters:*/
  double xi_low   = 1.0, dxi = 0.1; /*we will eventually want to control the range rather than the Bose gas parameter*/
  double fill     = 0.6;
  double fillingtolerance = 1e-4;
  double t2       = 0.0;
  int    N        = 100; 
  double Ndouble  = (double)N;
  int    Nxi      = 100; 
  
  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;
  

  /*function variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  
  double gapintegrand;
  double gapintegral;
  double D_maxk;

  double k;
  double kprime;
  int check = 0;

  /*vectors and matrices:*/
  gsl_vector *D             = gsl_vector_calloc(N);
  gsl_matrix *WFF0matrix    = gsl_matrix_calloc(N, N);

  /*for the chemical potential*/
  double mu_low = -3.0, mu, dmu = 0.001;
  int Nmu       = 6000;
  double muintegrand, muintegral = 0.0, fillingnum = 0.0;
  double mu_min_value;
  
  mu = -2.5 + 2.0 * 2.5 * fill;


  double mu_guess; 
  gsl_vector *mu_vector     = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu   = gsl_vector_calloc(Nmu);
  gsl_vector *fillingvector = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double)imu * dmu + mu_low;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }
  

  /*variables for ground state energy*/
  double E0integrand; 
  double E0 = 0.0;

  /* start of xi loop*/
  double xi, nB, rBB;
  int checkxi = 0; 

  

  printf("%s \t %s \t %s \n", "xi", "mu", "E0");

  for (int jxi = 0; jxi < Nxi; ++jxi)
  {
    check = 0;
    xi = xi_low + (double)jxi * dxi;
    
    /*derived parameters:*/
    nB       = constant * 2.0 * xi / (strength * lt); 
    rBB      = pow(lt / (2.0 * xi), 2.0); /*(nB * aBB^3)^(1/3) < 0.03 atmost! (1 percent depletion) */ 

    /*start guess for Delta and interaction:*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set( D, i, Deltaguess(k) );
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, rBB, rBF, nB, mB, N, eps0, lt) );
      }
    }

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
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
          Deltakprime   = gsl_vector_get(D, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          gapintegral  += gapintegrand; 
        }
        gsl_vector_set(D, i, gapintegral);
      }
      for (int imu = 0; imu < Nmu; ++imu)
      {
        mu_guess =  (double)imu * dmu + mu_low;
        muintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu_guess;
          Deltakprime   = gsl_vector_get(D, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          muintegrand   = 1.0 / (2.0 * Ndouble) * (1.0 -  epsilonkprime / EFkprime);
          muintegral   += muintegrand;
        }
        gsl_vector_set(find_min_mu, imu, pow(fill - muintegral, 2.0) );
        gsl_vector_set(fillingvector, imu, muintegral );
      }
      mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));
      fillingnum = gsl_vector_get(fillingvector, gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-4 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-4 && pow(fill - fillingnum, 2.0) < fillingtolerance )
      {
        break;
      }

      if ( fabs(mu) < 1e-4 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-4 && pow(fill - fillingnum, 2.0) < fillingtolerance )
      {
        break;
      }
      if (gsl_vector_max(D) < 1e-4 )
      {
        break;
      }
      mu = mu_min_value;
      ++check; 
    }

    /*Ground state energy: */
    E0 = 0.0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
      Deltakprime   = gsl_vector_get(D, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

      E0integrand   = - 1.0 / (4.0 * Ndouble) * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
      E0           += E0integrand; 
    }

    printf("%lg \t %lg \t %lg \n", xi, mu, E0 + fillingnum * mu);
    fprintf(stderr, "checkxi = %i, \t check = %i \n", ++checkxi, check );
  }
 
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  gsl_vector_free(fillingvector);
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  return 0;
}