#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_min.h>
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

struct T_and_mu   {double *T, *mu;};

double mufreeintegrand(double E, void *params) 
{
  struct T_and_mu *tmu = (struct T_and_mu*)params;
  
  double *T     = tmu->T;
  double *mu  = tmu->mu;
  double f = 1.0/2.0 * 1.0 / ( sqrt(E) * (exp( (E - *mu)/(*T) ) + 1.0) );
  return f;
}

double mufreeintegral(void *params) 
{
  gsl_function F = {.function = mufreeintegrand, .params = params};

  double epsabs = 1e-4, epsrel = 1e-4, result, err;

  int limit = 200;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);

  gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, workspace, &result, &err);
  gsl_integration_workspace_free(workspace);

  return result;
}

double master (double mu, void * T)
{
  struct T_and_mu tmu;

  tmu.T     = (double *)T;
  tmu.mu    = &mu;   

  return pow(1.0 - mufreeintegral( &tmu ) , 2.0);
}

int
main (void)
{
  /*variables */
  double rBB = 0.01;
  double rBF_high = 0.1, rBF_low = 0.05, drBF = 0.0001;    /*(nB * aBF^3)^(1/3) */
  double mB  = 7.0/40.0; /*mB/mF*/

  /*nB:*/
  /* double inverseBdist = 0.036; 
  double nB           = pow(inverseBdist, -3.0); */
  double nB  = 100.0; 
  double aB       = M_PI * rBB / pow(nB, 1.0/3.0);
  double retneg   = pow(1.0/mB, 2.0) * nB * aB;

  if (retneg < 20.0)
  {
    fprintf(stderr, "%s\n", "retneg too small" );
    return 0;
  }


  /*k-values:*/
  double k_up = 80.0, dk = 0.17;
  int N = (int) k_up/dk;
 
  /*L- matrix and eigenvalue, eigenvector matrices:*/
  gsl_matrix*                  L = gsl_matrix_calloc(N, N);
  gsl_eigen_nonsymm_workspace* w =  gsl_eigen_nonsymm_alloc (N);
  gsl_vector_complex*   eigenval = gsl_vector_complex_alloc(N);

  gsl_vector_view eigenval_real;
  gsl_vector_view eigenval_imag;

  /* Temperatures:*/
  double T_low, T_high, dT = 0.0001, T_guess = 0.13;
  printf("%s \t %s \t %s \n", "rBF", "largest eigenvalue", "TC");
  printf("\n \n");

  int check = 0, iter_number = 0;
  double TC = 0.0;
  double largest;
  double rBF;


  /*Free gas chemical potential*/
  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *M;
  gsl_min_fminimizer *s;
  
  double mufree_low = 0.9, mufree_up = 1.5, mufree_guess = 1.0;
  double mufree_low_2, mufree_up_2;
  double mu;
  gsl_function Fmin;

  Fmin.function = &master;

  M = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (M);


  for (rBF = rBF_high; rBF > rBF_low; rBF -= drBF)
  {
    T_low  = T_guess - 100.0 * dT;
    T_high = T_guess + dT; 

    check = 0;
    for (double T = T_high; T > T_low; T-=dT)
    {
      Fmin.params = &T;
      gsl_min_fminimizer_set (s, &Fmin, mufree_guess, mufree_low, mufree_up);
      iter = 0; 
      do
      {
        iter++;
        status            = gsl_min_fminimizer_iterate   (s);
        mufree_guess      = gsl_min_fminimizer_x_minimum (s);
        mufree_low_2      = gsl_min_fminimizer_x_lower   (s);
        mufree_up_2       = gsl_min_fminimizer_x_upper   (s);
        
        status            = gsl_min_test_interval (mufree_low_2, mufree_up_2, 0.001, 0.0);
        
        if (status == GSL_SUCCESS)
          mu = mufree_guess;
      }

    while (status == GSL_CONTINUE && iter < max_iter);

      largest = largesteigenval(L, dk, k_up, rBB, rBF, nB, mB, mu, T, eigenval, w, eigenval_real, eigenval_imag);
      if (largest > 1.0)
      {
        TC = T;
        break;
      }

      ++check;
    }

    fprintf(stderr, "largest = %lg, check = %i, iter = %i \n", largest, check, ++iter_number);

    if (TC < 1.0e-4 || check == 100)
    {
      break;
    }

    T_guess = TC;
    printf("%lg \t %lg \t %lg \n", rBF, largest, TC);
  }
  printf("%lg \t %lg \t %lg \n", rBF, largest, 0.0);

  


  gsl_vector_complex_free (eigenval);
  gsl_matrix_free(L);
  gsl_eigen_nonsymm_free (w);
  return 0;
}