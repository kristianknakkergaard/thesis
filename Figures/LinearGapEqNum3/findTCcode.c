#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_min.h>
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

gsl_matrix * Lmatrix(gsl_matrix *L, double dk, double k_up, double rBB, double rBF, double nB, double mB, double mu, double T)
{
  int N = (int) k_up/dk;
  double Lelement, k, kprime, epsilonkprime;

  double aB     = M_PI * rBB / pow(nB, 1.0/3.0);
  double retneg = pow(1.0/mB, 2.0) * nB * aB; 

  for (int i = 0; i < N; ++i)
  {
    if (retneg < 20)
    {
      fprintf(stderr, "%s\n", "retneg too small" );
      break;
    }
    k = ((double)i) * dk;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk;
      epsilonkprime = kprime * kprime - mu; /*mu set to 1!!! */
      Lelement      = -1.0/M_PI * WFF0(k, kprime, rBB, rBF, nB, mB) * tanh(epsilonkprime / (2.0 * T)) / (2.0 * epsilonkprime) * dk;
      gsl_matrix_set( L, i, iprime, Lelement );
    }
  }
  return L;
}

double largesteigenval(gsl_matrix *L, double dk, double k_up, double rBB, double rBF, double nB, double mB, double mu, double T, gsl_vector_complex* eigenval, gsl_eigen_nonsymm_workspace* w, gsl_vector_view eigenval_real, gsl_vector_view eigenval_imag)
{
  L = Lmatrix(L, dk, k_up, rBB, rBF, nB, mB, mu, T);

  gsl_eigen_nonsymm(L, eigenval, w);

  eigenval_real = gsl_vector_complex_real(eigenval);
  eigenval_imag = gsl_vector_complex_imag(eigenval);
  if (gsl_vector_max(&eigenval_imag.vector) > 0.1 )
  {
    return 10000.0;
  }
  return gsl_vector_max(&eigenval_real.vector);
}

struct myparams {gsl_matrix *L; double* dk, *k_up, *rBB, *rBF, *nB, *mB, *mu; gsl_vector_complex* eigenval; gsl_eigen_nonsymm_workspace* w; gsl_vector_view eigenval_real, eigenval_imag;};

double TCmaster (double T, void * params)
{
  struct myparams *myparams = (struct myparams*) params;

  gsl_matrix *L = myparams->L;
  double    *dk = myparams->dk;
  double  *k_up = myparams->k_up;
  double   *rBB = myparams->rBB;
  double   *rBF = myparams->rBF;
  double    *nB = myparams->nB;
  double    *mB = myparams->mB;
  double    *mu = myparams->mu;

  return fabs ( 1.0 - largesteigenval(L, *dk, *k_up, *rBB, *rBF, *nB, *mB, *mu, T, myparams->eigenval, myparams->w, myparams->eigenval_real, myparams->eigenval_imag) );
}

int
main (void)
{
  /*variables */
  double rBB = 0.01;
  double mB  = 7.0/40.0; /*mB/mF*/
  double mu  = 1.0;

  double rBF_high = 0.07, rBF_low = 0.049, drBF = 0.0001;    /*(nB * aBF^3)^(1/3) */
  

  /*nB:*/
  double inverseBdist = 0.036;
  double nB           = pow(inverseBdist, -3.0); 
  /* double nB  = 100.0; */


  /*k-values:*/
  double k_up = 60.0, dk = 0.27;
  int N = (int) k_up/dk;
 
  /*myparams */
  struct myparams myparams;
  myparams.dk   = &dk;
  myparams.mB   = &mB;
  myparams.k_up = &k_up;
  myparams.rBB  = &rBB;
  myparams.nB   = &nB;
  myparams.mu   = &mu;

  /* L-matrix and eigenvalue, eigenvector matrices:*/
  myparams.L        = gsl_matrix_calloc(N, N);
  myparams.w        = gsl_eigen_nonsymm_alloc (N);
  myparams.eigenval = gsl_vector_complex_alloc(N);

  gsl_vector_view eigenval_real;
  myparams.eigenval_real = eigenval_real;
  
  gsl_vector_view eigenval_imag;
  myparams.eigenval_imag = eigenval_imag;

  /*For minimization*/
  int status;
  int iter, max_iter = 100;
  const gsl_min_fminimizer_type *M;
  gsl_min_fminimizer *s;


  M = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (M);

  gsl_function Fmin;
  Fmin.function = &TCmaster;

  /* temperatures:*/
  double T_low = 0.001, T_up = 1.0, T_guess = 0.01295;
  double T_low_2, T_up_2;
  printf("%s \t %s \t %s \t %s \t %s \n", "rBF", "TC", "T_low", "T_up", "iter");
  printf("\n \n");

  int check = 0;
  double rBF;
  for (rBF = rBF_high; rBF > rBF_low; rBF -= drBF)
  {
    myparams.rBF = &rBF;

    Fmin.params = (void *)&myparams;
    gsl_min_fminimizer_set (s, &Fmin, T_guess, T_low, T_up);

    iter = 0;
    do
    {
      iter++;
      status       = gsl_min_fminimizer_iterate   (s);
      T_guess      = gsl_min_fminimizer_x_minimum (s);
      T_low_2      = gsl_min_fminimizer_x_lower   (s);
      T_up_2       = gsl_min_fminimizer_x_upper   (s);

      fprintf(stderr, "%lg \t %lg\n", T_low_2,T_up_2 );
      
      status 
        = gsl_min_test_interval (T_low_2, T_up_2, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        printf("%lg \t %lg \t %lg \t %lg \t %i \n", rBF, T_guess, T_low_2, T_up_2, iter);
    }

    while (status == GSL_CONTINUE && iter < max_iter);

    fprintf(stderr, "%i\n", ++check );
  }
  printf("%lg \t %lg \n", rBF, 0.0);

  


  gsl_vector_complex_free (myparams.eigenval);
  gsl_matrix_free(myparams.L);
  gsl_eigen_nonsymm_free (myparams.w);
  return status;
}