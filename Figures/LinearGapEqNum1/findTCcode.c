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

  double factor = 4.0 / (M_PI * M_PI) * (1.0/mB + mB + 2.0) * nB * pow(aBF, 2.0);
  double f      = - factor * log( ( pow(k+kprime, 2.0) + 2.0/(xi * xi) )/( pow(k-kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}


int
main (void)
{
  /*variables */
  double rBB = 0.01;
  double rBF = 0.106;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  double inverseBdist = 0.036;
  double nB           = pow(inverseBdist, -3.0); 
  /* double nB  = 100.0; */

  /*dependent variables:*/
  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);
  double retneg = pow(1.0/mB, 2.0) * nB * aB;

  /*k-values:*/
  double k_low = 0.0, k_up = 150.0, dk = 0.2;
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

  gsl_eigen_nonsymm_workspace* w =  gsl_eigen_nonsymm_alloc (N);
  gsl_vector_complex*   eigenval = gsl_vector_complex_alloc(N);

  gsl_vector_view eigenval_real;

  /* temperatures:*/
  double T_low = 0.17, T_high = 0.2, dT = 0.001;
  printf("%s \t %s \t %s \n", "T", "largest eigenvalue real", "largest eigenvalue imag");
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
        epsilonkprime = kprime * kprime - 0.649; /*mu set to 1!!! */
        Lelement      = -1.0/M_PI * gsl_matrix_get(W, i, iprime) * tanh(epsilonkprime / (2.0 * T)) / (2.0 * epsilonkprime) * dk;
        gsl_matrix_set( L, i, iprime, Lelement );
      }
    }
    gsl_eigen_nonsymm(L, eigenval, w); 

    eigenval_real = gsl_vector_complex_real(eigenval);

    printf( "%lg \t %lg  \n", T, gsl_vector_max(&eigenval_real.vector));

    if (gsl_vector_max(&eigenval_real.vector) < 1.0)
    {
      TC = T;
      break;
    }
    fprintf(stderr, "%i\n", ++check);
  }

  fprintf(stderr, "\n \n TC = %lg\n \n ", TC); 

  gsl_vector_complex_free (eigenval);
  gsl_matrix_free(L);
  gsl_eigen_nonsymm_free (w);
  return 0;
}