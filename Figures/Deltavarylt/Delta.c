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
double F(double q, double lt) /*q in units of kF*/
{
  double BBscatlength = 0.3; /*kF * aB*/
  double densityratio = 1.0;
  double xi           = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/
  return lt * lt/2.0 * ( q * q + 2.0/(xi * xi) );  
}

/*SET mF/mB, nB/nF^3 and kF * aBF */
struct q_and_lt {double *q, *lt;};

double Vintegrand0(double u, void *params)
{
  struct q_and_lt* qlt = ( struct q_and_lt* )params;
  double *q = qlt->q;
  double *lt = qlt->lt;

  double massratio    = 0.8; /*mF/mB*/
  double BFscatlength = 0.3; /*kF * aBF */
  double densityratio = 1.0; /*nB/nF^3*/
  double factor       = pow(2.0,5.0)/(M_PI*M_PI) * pow(BFscatlength, 2.0) * densityratio * (massratio + 1.0/massratio + 2.0);
  
  

  double integrand = -factor / (u + F(*q,*lt)) * exp(-u);
  return integrand;
}

double VFF0(void *params) 
{
  gsl_function Func = {.function = Vintegrand0, .params = params};

  double epsabs = 1e-6, epsrel = 1e-6, result, err;

  int limit = 200;
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(limit);
  gsl_integration_qagiu(&Func, 0.0, epsabs, epsrel, limit, workspace, &result, &err);
  
  gsl_integration_workspace_free(workspace);
  return result;
}

double WFF0(double k, double kprime, double lt)
{
  struct q_and_lt q_and_ltdiff;
  double diff = k - kprime;
  q_and_ltdiff.q   = &diff;
  q_and_ltdiff.lt  = &lt;

  struct q_and_lt q_and_ltsum;  
  double sum = k + kprime;
  q_and_ltsum.q  = &sum;
  q_and_ltsum.lt = &lt;
  

  double f = VFF0(&q_and_ltdiff) - VFF0(&q_and_ltsum);
  return f;
}

double Deltaguess(double k)
{
  return 0.4 * k / (pow(k,4) + 1);
}


int
main (void)
{
  /*k-values:*/
  double k_low = 0.0, k_up = 50.0, dk = 0.1;
  int N = (int) (k_up - k_low)/dk;
  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double integrand;
  double integral;
  double D_max;

  /*lt's*/
  double lt_low = 0.0003, lt_up = 0.3, dlt = 0.003;

  gsl_vector *D = gsl_vector_calloc(N);
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  /*For calculating mu:*/
  double mu_low = 0.97, dmu = 0.001;
  double muintegral, muintegrand;
  double mu_min_value = 0.98;
  int Nmu = 20;
  double mu, mu_guess;
  int mu_check;

  gsl_vector *mu_vector = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }


  double k;
  double kprime;
  /*start guess:*/
  for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      gsl_vector_set(D, i, Deltaguess(k));
    }

  /*Number of iterations*/
  int iter_number = 0;

  for (double lt = lt_low; lt < lt_up; lt+=dlt)
  {
    mu_check = 0;

    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, lt) );
      }
    }

    for (int iter = 0; iter < 50; ++iter)
    {
      mu    = mu_min_value;
      D_max = gsl_vector_max(D);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        integral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = kprime * kprime - mu;
          Deltakprime   = gsl_vector_get(D, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + fabs(Deltakprime * Deltakprime), 1.0 / 2.0);
          integrand     =  -1.0 / M_PI * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          integral     += integrand*dk; 
        }
        gsl_vector_set(D,i,integral);
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
      if ( fabs(mu - mu_min_value)/mu_min_value < 1e-2 && fabs(gsl_vector_max(D) - D_max)/D_max < 1e-3 )
      {
        break;
      }
      ++mu_check;
    }
    fprintf(stderr, "%i \t %i\n", ++iter_number, mu_check);
    printf("%lg \t %lg \t %lg \n", lt, gsl_vector_max(D), mu_min_value);
  }
  
  gsl_vector_free(D);
  gsl_matrix_free(WFF0matrix);
  gsl_vector_free(find_min_mu);
  gsl_vector_free(mu_vector);
  
  return 0;
}