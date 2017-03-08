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
  double xi     = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/ 

  double factor = 4.0 * (mB + 1.0/mB + 2.0) * pow(nB, 1.0 / 3.0) * rBF * rBF;
  double f      = - factor * log( ( pow( k + kprime, 2.0 ) + 2.0 / (xi * xi) ) / ( pow( k - kprime, 2.0 ) + 2.0 / (xi * xi) ) );
  return f; 
}

double Deltaasymp(double D_maxkT, double T, double TC)
{
  return D_maxkT * pow(1.0 - pow(T/TC, 3.0) , 1.0/2.0);
}

double Deltaguess(double k)
{
  return k / ( k * k + 3.3);
}



int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.1;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  double nB = 100.0;

  /*k-values:*/
  double k_low = 0.0, k_up = 50.0, dk = 0.01;
  int N = (int) (k_up - k_low)/dk;


  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double gapintegrand;
  double gapintegral;

  /*vectors and matrices:*/
  gsl_vector *D = gsl_vector_calloc(N);
  double D_maxk;
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;
  int check = 0;
  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k));

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k,kprime, rBB, rBF, nB, mB) );
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
        epsilonkprime = kprime * kprime - mu;
        Deltakprime   = gsl_vector_get(D, iprime);
        EFkprime      = sqrt(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime);
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
        EFkprime      = sqrt(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime);
        muintegrand   = 1.0/2.0 * (1.0 -  epsilonkprime / EFkprime);
        muintegral   += muintegrand*dk;
      }
      gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
    }
    mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

    if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
    {
      break;
    }
    mu = mu_min_value;
    ++check; 
  }

  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);
  double aBF    = M_PI * rBF / pow(nB, 1.0/3.0);
  double retneg = pow(1.0/mB, 2.0) * nB * aB;
  fprintf(stderr, "(nBaB^3)^(1/3) = %lg, (nBaBF^3)^(1/3) = %lg, mB/mF = %lg, nB/nF^3 = %lg, kF*aB = %lg, kF*aBF = %lg, (mF/mB)^2*nB/nF^3*kF*aB = %lg \n", rBB, rBF, mB, nB, aB, aBF, retneg);

  fprintf(stderr, "\n \n");
  fprintf(stderr, "%s \t %s \t %s \t %s \t %s\n", "x", "corfunc", "aveordparam", "T", "check" );
  fprintf(stderr, "\n \n");

  /*Real space correlation function*/

  /*x-values*/
  double x_up = 50.0, dx = 0.01; 
  double x0 = 2.0; 

  /*integral variables */
  double corfuncintegrand;
  double corfunc = 0.0;
  double aveordparamintegrand;
  double aveordparam = 0.0;

  for (double x = - x_up; x < x_up; x+=dx)
  {
    corfunc = 0.0;
    aveordparam = 0.0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime           = ((double)iprime) * dk + k_low;
      epsilonkprime    = kprime * kprime - mu;
      Deltakprime      = gsl_vector_get(D, iprime);
      EFkprime         = sqrt(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime);
      
      corfuncintegrand = 1.0 / M_PI * 1.0 / 2.0 * sin(kprime * x0) * sin(kprime * (x + x0)) * 1.0 / 2.0 * pow( 1.0 - epsilonkprime / EFkprime, 2.0);
      corfunc         += corfuncintegrand * dk;
      
      aveordparamintegrand = 1.0 / M_PI * sin(kprime * x) * Deltakprime / (2.0 * EFkprime);
      aveordparam         += aveordparamintegrand * dk;
    }
    fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %i \n", x, corfunc, aveordparam, 0.0, check);
  }
  fprintf(stderr, "\n \n" );
  
  /*T > 0*/

  const int numberT = 4; 
  double T[numberT] = {0.075, 0.100, 0.120, 0.127 };

  for (int iT = 0; iT < numberT; ++iT)
  {
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
          epsilonkprime = kprime * kprime - mu;
          Deltakprime   = gsl_vector_get(D, iprime);
          EFkprime      = sqrt(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime);
          gapintegrand  = -1.0/M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime) * tanh(EFkprime / (2.0 * T[iT]));
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
          EFkprime      = sqrt(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime);
          muintegrand   = epsilonkprime / EFkprime * ( 1.0/(exp(EFkprime / T[iT]) + 1.0) - 1.0/2.0) + 1.0/2.0;
          muintegral   += muintegrand*dk;
        }
        gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
      {
        break;
      }
      mu = mu_min_value;
      ++check; 
    }

    /*Real space correlation function*/

    for (double x = - x_up; x < x_up; x+=dx)
    {
      corfunc = 0.0;
      aveordparam = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime           = ((double)iprime) * dk + k_low;
        epsilonkprime    = kprime * kprime - mu;
        Deltakprime      = gsl_vector_get(D, iprime);
        EFkprime         = sqrt(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime);
        
        corfuncintegrand = 1.0 / M_PI * 1.0 / 2.0 * sin(kprime * x0) * sin(kprime * (x + x0)) * ( (1.0 + pow(epsilonkprime / EFkprime, 2.0) ) * pow( 1.0 / ( exp(EFkprime / T[iT]) + 1.0 ), 2.0) + 1.0 / 2.0 * pow( 1.0 - epsilonkprime / EFkprime, 2.0 ) * tanh(EFkprime / (2.0 * T[iT]))  );
        corfunc         += corfuncintegrand * dk;

        aveordparamintegrand =  1.0 / M_PI * sin(kprime * x) * Deltakprime / (2.0 * EFkprime) * tanh(EFkprime / (2.0 * T[iT]));
        aveordparam         += aveordparamintegrand * dk;

      }
      fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %i \n", x, corfunc, aveordparam, T[iT], check);
    }

    fprintf(stderr, "\n \n" );
  }
  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}