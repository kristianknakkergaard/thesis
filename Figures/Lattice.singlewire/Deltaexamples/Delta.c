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
  double xi     = lt / (2.0 * sqrt( rBB ) ); /*xi / a*/
  double factor = sqrt(2.0) * (mB + 1.0 / mB + 2.0) * rBF * rBF * 1.0 / (lt * lt * lt) * eps0 * 1.0 / nB * pow(rBB, -1.0 / 2.0);

  return - factor * exp(- sqrt(2.0) * fabs(x) / xi ); 
}

double WFF0(double k, double kprime, double rBB, double rBF, double nB, double mB, int N, double eps0, double lt)
{
  double W = 0.0;

  for (int i = 1; i < N; ++i)
  {
    double x  = (double) i;
    W += 1.0 / 2.0 * ( cos((kprime - k) * i ) - cos( (kprime + k) * i ) ) * Vindx( x, rBB, rBF, nB, mB, eps0, lt); 
  }

  return W; 
}

double Deltaguess(double k, double a, double b)
{
  return a * sin(k) + b * sin(2.0 * k);
}


int
main (void)
{
  /*restricted variables */
  double t1       = 1.0;  /*the unit of energy */
  double strength = 4.0; /*total strength of interaction*/
  double mB       = 7.0/40.0; /*mB/mF*/
  double rBF      = 0.005;    /*(nB * aBF^3)^(1/3) */
  double eps0     = 0.5; /*energy of particle box */
  double lt       = 0.1;
  double constant = sqrt(2.0) * (mB + 1.0 / mB + 2.0) * rBF * rBF * 1.0 / (lt * lt * lt) * eps0;  

  /*free parameters:*/
  double xi       = 5.0; /*we will eventually want to control the range rather than the Bose gas parameter*/
  double t2       = 1.0;
  int    N        = 200; 
  double Ndouble  = (double)N;
  double filling[3]  = {0.2, 0.5, 0.8};
  double filling_value;
  double filltolerance = 1e-4; 
  double fillingnum; 

  /*derived parameters:*/
  double rBB      = pow(lt / (2.0 * xi), 2.0); /*(nB * aBB^3)^(1/3) < 0.03 atmost! (1 percent depletion) */
  double nB       = constant * 1.0 / (strength * sqrt(rBB));


  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  
  double gapintegrand;
  double gapintegral;
  double D_maxk;
  double Delta1ki, Delta2ki;

  double k;
  double kprime;
  int check = 0; 

  /*vectors and matrices:*/
  gsl_vector *D1k          = gsl_vector_calloc(N);
  gsl_vector *D2k          = gsl_vector_calloc(N);
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);


  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, rBB, rBF, nB, mB, N, eps0, lt) );
    }
  }

  /*for the chemical potential*/

  double mu_low     = -3.0, dmu = 0.001;
  double mu         = 0.3;
  int Nmu           = 6000;
  double muintegral = 0.0, muintegrand;
  double mu_min_value;
  
  double mu_guess;
  gsl_vector *mu_vector     = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu   = gsl_vector_calloc(Nmu);
  gsl_vector *fillingvector = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }


  /*T = 0: */
  printf("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n", "k", "Delta1k", "Delta2k", "epsilonk", "site", "Delta1ki", "Delta2ki", "check");
  for (int i_filling = 0; i_filling <= 2; ++i_filling)
  {
    filling_value = filling[i_filling];
    /*start guess for Delta:*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set( D1k, i, Deltaguess(k, 1.0, 0.0) );
    }
    for (int j = 0; j < 1000; ++j)
    {
      D_maxk = gsl_vector_max(D1k);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
          Deltakprime   = gsl_vector_get(D1k, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          gapintegral  += gapintegrand; 
        }
        gsl_vector_set(D1k, i, gapintegral);
      }
      for (int imu = 0; imu < Nmu; ++imu)
      {
        mu_guess    = (double) mu_low + imu * dmu;
        muintegral  = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu_guess;
          Deltakprime   = gsl_vector_get(D1k, iprime );
          EFkprime      = pow( epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0 );
          muintegrand   = 1.0 / ( 2.0 * Ndouble ) * ( 1.0 -  epsilonkprime / EFkprime );
          muintegral   += muintegrand;
        }
        gsl_vector_set(find_min_mu, imu, pow(filling_value - muintegral, 2.0) );
        gsl_vector_set(fillingvector, imu, muintegral );
      }
      mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));
      fillingnum   = gsl_vector_get(fillingvector, gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value) / mu_min_value < 1e-4     && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D1k) - D_maxk)/D_maxk < 1e-4 )
      {
        break;
      }
      if ( fabs(mu - mu_min_value) < 1e-4 && fabs(mu) < 1e-4 && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D1k) - D_maxk)/D_maxk < 1e-4 )
      {
        break;
      }
      if ( fabs(mu - mu_min_value) / mu_min_value < 1e-4     && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D1k)) < 1e-4 )
      {
        break;
      }
      if ( fabs(mu - mu_min_value) < 1e-4 && fabs(mu) < 1e-4 && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D1k)) < 1e-4 )
      {
        break;
      }
      mu          = mu_min_value;
      ++check; 
    }

    /*start guess for Delta:*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set(D2k, i, Deltaguess(k, 0.0, 1.0) );
    }

    for (int j = 0; j < 1000; ++j)
    {
      D_maxk = gsl_vector_max(D2k);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
          Deltakprime   = gsl_vector_get(D2k, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          gapintegral  += gapintegrand; 
        }
        gsl_vector_set(D2k, i, gapintegral);
      }
      for (int imu = 0; imu < Nmu; ++imu)
      {
        mu_guess    = (double) mu_low + imu * dmu;
        muintegral  = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu_guess;
          Deltakprime   = gsl_vector_get(D2k, iprime );
          EFkprime      = pow( epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0 );
          muintegrand   = 1.0 / ( 2.0 * Ndouble ) * ( 1.0 -  epsilonkprime / EFkprime );
          muintegral   += muintegrand;
        }
        gsl_vector_set(find_min_mu, imu, pow(filling_value - muintegral, 2.0) );
        gsl_vector_set(fillingvector, imu, muintegral );
      }
      mu_min_value = gsl_vector_get(mu_vector, gsl_vector_min_index(find_min_mu));
      fillingnum   = gsl_vector_get(fillingvector, gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value) / mu_min_value < 1e-4     && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D2k) - D_maxk)/D_maxk < 1e-4 )
      {
        break;
      }
      if ( fabs(mu - mu_min_value) < 1e-4 && fabs(mu) < 1e-4 && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D2k) - D_maxk)/D_maxk < 1e-4 )
      {
        break;
      }
      if ( fabs(mu - mu_min_value) / mu_min_value < 1e-4     && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D2k)) < 1e-4 )
      {
        break;
      }
      if ( fabs(mu - mu_min_value) < 1e-4 && fabs(mu) < 1e-4 && pow(fillingnum - filling_value, 2.0) < filltolerance && fabs(gsl_vector_max(D2k)) < 1e-4 )
      {
        break;
      }
      mu          = mu_min_value;
      ++check; 
    }


    /*We list the function values:*/ 
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      epsilonkprime = - t1 * cos(k) - t2 * cos(2.0 * k) - mu;
      Delta1ki = 0.0, Delta2ki = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        Delta1ki += 1.0 / (Ndouble) * (gsl_vector_get(D1k, iprime) * sin(kprime * i));
        Delta2ki += 1.0 / (Ndouble) * (gsl_vector_get(D2k, iprime) * sin(kprime * i));
      }
      printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %i \n", k, gsl_vector_get(D1k, i), gsl_vector_get(D2k, i), epsilonkprime, (double)i, Delta1ki, Delta2ki, check);
    }
    printf("\n \n");
  }
  
  
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  gsl_vector_free(fillingvector);
  gsl_vector_free(D1k);
  gsl_vector_free(D2k);
  gsl_matrix_free(WFF0matrix);
  return 0;
}