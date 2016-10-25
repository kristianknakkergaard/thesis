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
double WFF0(double k, double kprime)
{
  double BBscatlength = 0.3; /*kF * aB*/
  double densityratio = 1.0;   /*nB/nF^3 */
  double xi           = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/ 
  double massratio    = 0.8; /*mB/mF*/
  double BFscatlength = 0.3; /*kF * aBF */
  double factor       = pow(2.0,5.0)/(M_PI*M_PI) * pow(BFscatlength, 2.0) * densityratio * (massratio + 1.0/massratio + 2.0);

  double f = - factor * log( ( pow(k+kprime, 2.0) + 2.0/(xi * xi) )/( pow(k-kprime, 2.0) + 2.0/(xi * xi) ) );
  return f;
}


double Deltaguess(double k)
{
  return 0.4 * k / (pow(k,4) + 1);
}

double Deltaasymp(double D_maxkT, double T, double TC)
{
  return D_maxkT * pow(1.0 - pow(T/TC, 3.0) , 1.0/2.0);
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

  double epsabs = 1e-6, epsrel = 1e-6, result, err;

  int limit = 400;
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
  /*k-values:*/
  double k_low = 0.0, k_up = 50.0, dk = 0.01;
  int N = (int) (k_up - k_low)/dk;
  /*fermion mass*/

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double gapintegrand;
  double gapintegral;
  double mu = 1.0;

  gsl_vector *D = gsl_vector_calloc(N);
  double D_maxk;
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  double k;
  double kprime;

  int check = 0;

  /*For calculating mu:*/
  double mu_low = 0.982, dmu = 0.0001;
  double muintegral, muintegrand;
  double mu_min_value;
  int Nmu = 350;
  double mu_guess;
  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  /*start guess for Delta:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D, i, Deltaguess(k));

    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime = ((double)iprime) * dk + k_low;
      gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k,kprime) );
    }
  }
  
  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }
  
  /*T = 0: */
  for (int j = 0; j < 50; ++j)
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
      gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
    }
    mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

    if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
    {
      break;
    }
    mu = mu_min_value;
    ++check; 
  }

  double D_maxkT = D_maxk; 

  double kmax = (double)gsl_vector_max_index(D) * dk + k_low;
  fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %lg \t %i \n", 0.0, kmax, gsl_vector_max(D), mu, 1.0, check);


  /*We list the function values:*/
  for (int i = 0; i < N-1; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), 0.0);
  }
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    printf("%lg \t %lg \t %lg \n", k, gsl_vector_get(D,i), 0.0);
  }
  printf("\n\n");

  /*§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§*/
  /*Now: temperature T>0 */
  double T_low = 0.001, T_high = 0.2, dT = 0.001;
  int kmax_index;

  /*§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§*/
  /*Free gas chemical potential*/

  int status;
  int iter = 0, max_iter = 100;
  const gsl_min_fminimizer_type *L;
  gsl_min_fminimizer *s;
  
  double mufree_low = 0.9, mufree_up = 1.5, mufree_guess = 1.0;
  double mufree_low_2, mufree_up_2;
  gsl_function Fmin;

  Fmin.function = &master;

  L = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc (L);

  for (double T = T_low; T < T_high; T+=dT)
  {
    check = 0;
    for (int j = 0; j < 200; ++j)
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
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  =  -1.0/M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0*EFkprime) * tanh(EFkprime / (2.0*T));
          gapintegral  += gapintegrand*dk; 
        }
        gsl_vector_set(D, i, gapintegral);
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
          muintegrand   = epsilonkprime / EFkprime * ( 1.0/(exp(EFkprime / T) + 1.0) - 1.0/2.0) + 1.0/2.0;
          muintegral   += muintegrand*dk;
        }
        gsl_vector_set(find_min_mu, imu, fabs(1.0 - muintegral) );
      }
      mu_min_value = gsl_vector_get(mu_vector,gsl_vector_min_index(find_min_mu));

      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-3 && fabs(gsl_vector_max(D) - D_maxk)/D_maxk < 1e-3 )
      {
        break;
      }
      mu = mu_min_value;
      ++check;
    }

    if (gsl_vector_max(D) < 1e-3)
    {
      break;
    }
    kmax_index  = gsl_vector_max_index(D);
    kmax        = ((double)kmax_index) * dk + k_low;

    Fmin.params = &T;
    gsl_min_fminimizer_set (s, &Fmin, mufree_guess, mufree_low, mufree_up);

    iter = 0; 
    do
    {
      iter++;
      status            = gsl_min_fminimizer_iterate (s);
      mufree_guess      = gsl_min_fminimizer_x_minimum (s);
      mufree_low_2      = gsl_min_fminimizer_x_lower (s);
      mufree_up_2       = gsl_min_fminimizer_x_upper (s);
      
      status 
        = gsl_min_test_interval (mufree_low_2, mufree_up_2, 0.001, 0.0);

      if (status == GSL_SUCCESS)
        fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %lg \t %i \n", T, kmax, gsl_vector_max(D), mu_min_value, mufree_guess, check);
    }

    while (status == GSL_CONTINUE && iter < max_iter);

    for (int i = 0; i < N-1; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg \n", k-k_up+dk, -gsl_vector_get(D,N-1-i), T, mu_min_value);
    }
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      printf("%lg \t %lg \t %lg \t %lg  \n", k, gsl_vector_get(D,i), T, mu_min_value );
    }

    printf("\n\n");
  }

  double TC = 0.1504; 
  fprintf(stderr, "\n\n");

  for (double T = 0; T < T_high; T+=dT)
  {
    fprintf(stderr, "%lg \t %lg \n", T, Deltaasymp(D_maxkT, T, TC) ); 
  }

  gsl_min_fminimizer_free(s);
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  return status;
}