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

double Deltaguess(double k)
{
  return 0.4 * k / (pow(k,4) + 1);
}


int
main (void)
{
  /*variables */
  double rBB = 0.01; /*(nB * aBB^3)^(1/3) <= 0.03 atmost! (1 percent depletion) */
  double rBF = 0.04;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/
  double nB  = 100.0;  

  /*k-values:*/
  double k_low = 0.0, k_up = 50.0, dk = 0.01;
  int N = (int) (k_up - k_low)/dk;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double integrand;
  double integral;

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

  /*Temperatures:*/
  double T_low = 0.0, T_high = 0.5, dT = 0.0001;

  /*k for maximum of Delta*/
  int kmax_index;
  double kmax;

  for (double T = T_low; T < T_high; T+=dT)
  {
    check = 0;
    for (int iter = 0; iter < 300; ++iter)
    {
      D_maxk = gsl_vector_max(D);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        integral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = kprime * kprime - 1.0;
          Deltakprime   = gsl_vector_get(D,iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + fabs(Deltakprime * Deltakprime), 1.0/2.0);

          if (T == 0.0)
          {
            integrand   =  -1.0/M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0*EFkprime);
          }
          else
            integrand   =  -1.0/M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0*EFkprime) * tanh(EFkprime / (2.0*T));

          integral  += integrand*dk;
        }
        gsl_vector_set(D, i, integral);
      }
      if ( fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
      {
        break;
      }
      ++check;
    }
    if (gsl_vector_max(D) < 1e-2)
    {
      break;
    }
    kmax_index  = gsl_vector_max_index(D);
    kmax        = ((double)kmax_index) * dk + k_low;
    /*T-dependency:*/
    fprintf(stderr, "%lg \t %lg \t %lg \t %i \n", T, kmax, gsl_vector_max(D), check);

    /*k-dependency:*/
    for (int i = 0; i < N-1; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg  \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), T);
    }
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg  \n", k, gsl_vector_get(D,i), T );
    }
    printf("\n\n");
  }
  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  return 0;
}