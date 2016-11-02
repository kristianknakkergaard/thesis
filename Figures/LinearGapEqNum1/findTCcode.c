#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
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

  double factor = 4.0 / (M_PI*M_PI) * (1.0/mB + 1.0) * nB * pow(aBF, 2.0);
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
  double rBB = 0.01;
  double rBF = 0.075;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  double inverseBdist = 0.036;  
  double nB  = pow(inverseBdist, -3.0); 

  /*dependent variables:*/
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);
  double retneg = pow(1.0/mB, 2.0) * nB * aB;

  /*k-values:*/
  double k_low = 0.0, k_up = 70.0, dk = 0.145;
  int N = (int) (k_up - k_low)/dk;

  /*to calculate gap:*/
  double epsilonkprime;
 
  /*L- matrix and eigenvalue, eigenvector matrices:*/
  gsl_matrix *L = gsl_matrix_calloc(N, N);
  double Lelement;

  gsl_matrix *W = gsl_matrix_calloc(N, N);
  double k;
  double kprime;
  for (int i = 0; i < N; ++i)
    {
      if (retneg < 20)
      {
        fprintf(stderr, "%s\n", "retneg too small" );
        break;
      }

      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        gsl_matrix_set( W, i, iprime, WFF0(k,kprime, rBB, rBF, nB, mB) );
      }
    }

  gsl_eigen_symm_workspace* w =  gsl_eigen_symm_alloc (N);
  gsl_vector* eigenval        = gsl_vector_alloc(N);

  /* temperatures:*/
  double T_low = 0.05, T_high = 3.0, dT = 0.01;
  printf("%s \t %s \n", "T", "largest eigenvalue");
  printf("\n \n");

  int check = 0;
  double TC = 0.0;
  for (double T = T_low; T < T_high; T+=dT)
  {
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = kprime * kprime - 1.0; /*mu set to 1!!! */
        Lelement      = -1.0/M_PI * gsl_matrix_get(W, i, iprime) * tanh(fabs(epsilonkprime) / (2.0 * T))/(2.0 * fabs(epsilonkprime)) * dk;
        gsl_matrix_set( L, i, iprime, Lelement );
      }
    }
    gsl_eigen_symm(L, eigenval, w);
    printf("%lg \t %lg \t %lg \n", T, gsl_vector_max(eigenval), gsl_vector_min(eigenval));

    if (gsl_vector_max(eigenval) < 1.0)
    {
      TC = T;
      break;
    }
    fprintf(stderr, "%i\n", ++check);
  }

  fprintf(stderr, "\n \n TC = %lg\n \n ", TC); 

  gsl_vector_free (eigenval);
  gsl_matrix_free(L);
  gsl_eigen_symm_free (w);
  return 0;
  fprintf(stderr, "iter_number = %i\n", ++check);
}