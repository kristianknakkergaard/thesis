#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
/*Units: */
/*Energy: the Fermi energy \epsilon_F */
/*Momentum: the fermi momentum k_F */
/*Length: */
/*All is for T=0*/

/*SET lt, nB/nF^3, kF * aB HERE*/
double F(double q) /*q in units of kF*/
{
  double BBscatlength = 0.3; /*kF * aB*/
  double densityratio = 1;
  double xi           = M_PI/sqrt(8.0 * densityratio * BBscatlength ); /*xi * kF*/ 
  double lt           = 0.5 * xi; /*lt * kF*/
  return lt * lt/2.0 * ( q * q + 2.0/(xi * xi) );  
}

/*SET mF/mB, nB/nF^3 and kF * aBF */
double Vintegrand0(double u, void *params) 
{
  double massratio    = 0.8; /*mF/mB*/
  double BFscatlength = 0.3; /*kF * aBF */
  double densityratio = 0.9; /*nB/nF^3*/
  double factor       = pow(2.0,5.0)/(M_PI*M_PI) * pow(BFscatlength, 2.0) * densityratio * (massratio + 1.0/massratio + 2.0);
  double q = *(double*)params; /* */
  double integrand = -factor / (u + F(q)) * exp(-u);
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

double WFF0(double k, double kprime)
{
  double sum = k + kprime;
  double diff = k - kprime;  

  double f = VFF0(&diff) - VFF0(&sum);
  return f;
}

double Deltaguess(double k)
{
  return 0.04 * k / (pow(k,4) + 1);
}


int
main (void)
{
  /*k-values:*/
  double k_low = 0.0, k_up = 25.0, dk = 0.05;
  int N = (int) (k_up - k_low)/dk;
  /*fermion mass*/

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  double integrand;
  double integral;

  gsl_vector *D = gsl_vector_calloc(N);


  double k;
  /*start guess:*/
  for (int i = 0; i < N; ++i)
  {
    k = ((double)i) * dk + k_low;
    gsl_vector_set(D,i,Deltaguess(k));
  }
  double kprime;
  for (int iter = 0; iter < 4; ++iter)
  {
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      integral = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        epsilonkprime = kprime * kprime - 1.0;
        Deltakprime = gsl_vector_get(D,iprime);
        EFkprime = pow(epsilonkprime * epsilonkprime + fabs(Deltakprime * Deltakprime), 1.0/2.0);
        integrand =  -1.0/M_PI * WFF0(k,kprime) * Deltakprime/(2.0 *EFkprime);
        integral+= integrand*dk; 
      }
      gsl_vector_set(D,i,integral);
    }
  }

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

  /*Now: temperature T>0 */
  double T_low = 0.001, T_high = 0.01, dT = 0.0001;
  int kmax_index;
  double kmax;

  for (double T = T_low; T < T_high; T+=dT)
  {
    for (int iter = 0; iter < 3; ++iter)
    {
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        integral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime = ((double)iprime) * dk + k_low;
          epsilonkprime = kprime * kprime - 1.0;
          Deltakprime = gsl_vector_get(D,iprime);
          EFkprime = pow(epsilonkprime * epsilonkprime + fabs(Deltakprime * Deltakprime), 1.0/2.0);
          integrand =  -1.0/M_PI * WFF0(k,kprime) * Deltakprime / (2.0*EFkprime) * tanh(EFkprime / (2.0*T));
          integral+= integrand*dk;
        }
        gsl_vector_set(D, i, integral);
      }
    }
    kmax_index  = gsl_vector_max_index(D);
    kmax        = ((double)kmax_index) * dk + k_low;

    fprintf(stderr, "%lg \t %lg \t %lg \n", T, kmax, gsl_vector_max(D));
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
  return 0;
}