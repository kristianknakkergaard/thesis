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

  double factor = pow(2.0,5.0)/(M_PI*M_PI) * pow(aBF, 2.0) * nB * (mB + 1.0/mB + 2.0);
  double f      = - factor * log( ( pow(k+kprime, 2.0) + 2.0/(xi * xi) )/( pow(k-kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}

double Deltaasymp(double D_maxkT, double T, double TC)
{
  return D_maxkT * pow(1.0 - pow(T/TC, 3.0) , 1.0/2.0);
}

double Deltaguess(double k)
{
  return 0.4 * k / (pow(k,4) + 1);
}


int
main (void)
{
  /*variables */
  double nB_low = 1e2, nB_up = 1e7, dnB = 1e4; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.04;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/
  double rBB = 0.01;

  /*dependent variables:*/
  double aB;
  double aBF;
  double retneg;

  /*k-values:*/
  double k_low = 0.0, k_up = 100.0, dk = 0.1;
  int N = (int) (k_up - k_low)/dk;

  /*to calculate gap:*/
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

  /*For calculating mu:*/
  double mu_low = 0.8, dmu = 0.0005;
  double muintegral, muintegrand;
  double mu_min_value = 1.0;
  double mu;
  int Nmu = 400;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k));
  }

  /*Number of iterations*/
  int iter_number = 0;
  
  printf("%s \t %s \t %s \t %s \t %s \t %s \n", "1/nB^(1/3)", "Delta", "mu", "aB", "retneg", "aBF");
  printf("\n \n");

  /*nB dependency for T = 0*/
  for (double nB = nB_up; nB > nB_low; nB-=dnB)
  {
    aB     = M_PI * rBB / pow(nB, 1.0/3.0);
    aBF    = M_PI * rBF / pow(nB, 1.0/3.0);
    retneg = pow(1.0/mB, 2.0) * nB * aB;

    check  = 0;
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k,kprime, rBB, rBF, nB, mB) );
      }
    }

    /*Find gap:*/
    for (int j = 0; j < 300; ++j)
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
        gsl_vector_set(find_min_mu, imu, pow(1.0 - muintegral, 2.0) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));
      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
      {
        break;
      }
      ++check; 
    }
    fprintf(stderr, "%i \t %i\n", ++iter_number, check);
    if (retneg < 20.0)
    {
      break;
    }
    printf("%.4lg \t %.4lg \t %.4lg \t %.4lg \t %.4lg \t %.4lg \n", pow(nB, -1.0/3.0), gsl_vector_max(D), mu_min_value, aB, retneg, aBF);
  }

  /*independen variables to print:*/
  
  fprintf(stderr, "(nBaBF^3)^(1/3) = %lg, mB/mF = %lg, (nBaB^3)^(1/3) = %lg \n", rBF, mB, rBB);


  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return 0;
}