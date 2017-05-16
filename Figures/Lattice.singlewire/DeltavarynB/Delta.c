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

double Vindx (double x, double rBB, double rBF, double nB, double mB, double eps0, double lt)
{
  double xi     = 1.0 / 2.0 * sqrt(1.0 / rBB ) * pow(nB, -1.0 / 3.0) ; /*xi * kF*/
  double factor = 2.0 * sqrt(2.0) * (mB + 1.0 / mB + 2.0) * rBF * rBF * pow(nB, 1.0 / 3.0) * xi * 1.0 / (lt * lt) * eps0;

  return - factor * exp(- sqrt(2.0) * fabs(x) / xi ); 
}

double WFF0(double k, double kprime, double rBB, double rBF, double nB, double mB, int N, double eps0, double lt)
{
  double W = 0.0;

  for (int i = 1; i <= N; ++i)
  {
    double x  = (double) i;
    W += 1.0 / 2.0 * ( cos((kprime - k) * i ) - cos( (kprime + k) * i ) ) * Vindx( x, rBB, rBF, nB, mB, eps0, lt); 
  }

  return W; 
}

double Deltaguess(double k, double ak, double b2k)
{
  return ak * sin(k) - b2k * sin(2.0 * k);
}


int
main (void)
{
  /*variables */
  double rBB     = 0.0001; /*(nB * aBB^3)^(1/3) < 0.03 atmost! (1 percent depletion) */
  double rBF     = 0.15;    /*(nB * aBF^3)^(1/3) */
  double nB_low  = 0.000001, nB_up = 0.0001, dnB = 0.000002; 
  double mB      = 7.0/40.0; /*mB/mF*/
  double lt      = 0.3;
  
  double t1      = 1.0;  /*the unit of energy, also gives unit of length */
  double t2      = 1.0;
  double eps0    = 0.02; 
  
  int    N       = 400; 
  double Ndouble = (double)N;

  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;

  /*variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  
  double gapintegrand;
  double gapintegral;
  double D_maxk;

  double k;
  double kprime;
  int check1k = 0, check2k = 0;

  /*vectors and matrices:*/
  gsl_vector *D1k        = gsl_vector_calloc(N);
  gsl_vector *D2k        = gsl_vector_calloc(N);
  gsl_matrix *WFF0matrix = gsl_matrix_calloc(N, N);

  double mu = 0.6;
  double muintegral1k = 0.0, muintegral2k = 0.0, muintegrand;

  /*variables for CS1 */
  double CS11k = 0;
  double CS12k = 0;
  double CS1integrand, dDelta;
  double depsilon; 

  /*variables for ground state energyv*/
  double E0integrand; 
  double E0integral1k = 0; 
  double E0integral2k = 0; 

  /*other parameters*/
  double xi; 
  double vFoverc0; 

  printf("%s \t %s \t %s \t %s \n", "k", "Delta1k", "Delta2k", "epsilonk");

  
  for (double nB = nB_low; nB < nB_up; nB+=dnB)
  {
    /*interaction and start guess for Deltas:*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      /*sin(k) only:*/
      gsl_vector_set( D1k, i, Deltaguess(k, 1.0, 0.0) );
      /*sin(2k) only:*/
      gsl_vector_set( D2k, i, Deltaguess(k, 0.0, 1.0) );

      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, rBB, rBF, nB, mB, N, eps0, lt) );
      }
    }
    
    /*sin(k) only: */
    check1k = 0; 
    for (int j = 0; j < 1000; ++j)
    {
      D_maxk = gsl_vector_max(D1k);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
          Deltakprime   = gsl_vector_get(D1k, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          gapintegral  += gapintegrand; 
        }
        gsl_vector_set(D1k, i, gapintegral);
      }
    
      if ( fabs(gsl_vector_max(D1k) - D_maxk)/D_maxk < 1e-4 )
      {
        break;
      }

      if (gsl_vector_max(D1k) < 1e-4 )
      {
        break;
      }
      ++check1k; 
    }

    /*CS-invariant:*/
    CS11k = 0; 
    muintegral1k = 0; 
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
      Deltakprime   = gsl_vector_get(D1k, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
      muintegrand   = 1.0 / (2.0 * Ndouble) * (1.0 -  epsilonkprime / EFkprime);
      muintegral1k += muintegrand;
      
      if (iprime < N - 1)
        {
          dDelta       = (gsl_vector_get(D1k, iprime + 1) - Deltakprime ) / dk;
          depsilon     = t1 * sin(kprime) + 2.0 * t2 * sin(2.0 * kprime);
          CS1integrand = 1.0 / (2.0 * M_PI) * (Deltakprime * depsilon - epsilonkprime * dDelta)/(EFkprime * EFkprime);
          CS11k       += dk * CS1integrand; 
        }
    }

    /*Ground state energy: */
    E0integral1k = 0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
      Deltakprime   = gsl_vector_get(D1k, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

      E0integrand   = - 1.0 / 4.0 * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
      E0integral1k += E0integrand*dk; 
    }

    /*sin(2k) only: */
    check2k = 0; 
    for (int j = 0; j < 1000; ++j)
    {
      D_maxk = gsl_vector_max(D2k);
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gapintegral = 0.0;
        for (int iprime = 0; iprime < N; ++iprime)
        {
          kprime        = ((double)iprime) * dk + k_low;
          epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
          Deltakprime   = gsl_vector_get(D2k, iprime);
          EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
          gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
          gapintegral  += gapintegrand; 
        }
        gsl_vector_set(D2k, i, gapintegral);
      }
    
      if ( fabs(gsl_vector_max(D2k) - D_maxk)/D_maxk < 1e-4 )
      {
        break;
      }

      if (gsl_vector_max(D2k) < 1e-4 )
      {
        break;
      }
      ++check2k; 
    }

    /*CS-invariant:*/
    CS12k         = 0;
    muintegral2k  = 0; 
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
      Deltakprime   = gsl_vector_get(D2k, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
      muintegrand   = 1.0 / (2.0 * Ndouble) * (1.0 -  epsilonkprime / EFkprime);
      muintegral2k   += muintegrand;
      
      if (iprime < N - 1)
        {
          dDelta       = (gsl_vector_get(D2k, iprime + 1) - Deltakprime ) / dk;
          depsilon     = t1 * sin(kprime) + 2.0 * t2 * sin(2.0 * kprime);
          CS1integrand = 1.0 / (2.0 * M_PI) * (Deltakprime * depsilon - epsilonkprime * dDelta)/(EFkprime * EFkprime);
          CS12k       += dk * CS1integrand; 
        }
    }

    /*Ground state energy: */
    E0integral2k = 0;
    for (int iprime = 0; iprime < N; ++iprime)
    {
      kprime        = ((double)iprime) * dk + k_low;
      epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu;
      Deltakprime   = gsl_vector_get(D2k, iprime);
      EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

      E0integrand   = - 1.0 / 4.0 * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
      E0integral2k += E0integrand*dk; 
    }
    

    xi       = sqrt(M_PI / ( 8.0 * rBB) ) * pow(nB, -1.0 / 3.0); 
    vFoverc0 = sqrt(M_PI) / 2.0 * mB * 1.0 / sqrt(rBB) * pow(nB, -1.0 / 3.0); 

    fprintf(stderr, "%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %i \t %i \n", nB, E0integral1k + mu, E0integral2k + mu, xi, vFoverc0, CS11k, CS12k, muintegral1k, muintegral2k, check1k, check2k);


    /*We list the function values:*/ 
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;
      epsilonkprime = - t1 * cos(k) - t2 * cos(2.0 * k) - mu;
      printf("%lg \t %lg \t %lg \t %lg \n", k, gsl_vector_get(D1k, i), gsl_vector_get(D2k, i), epsilonkprime);
    }
    printf("\n \n");
  }
  
  
  gsl_vector_free(D1k);
  gsl_vector_free(D2k);
  gsl_matrix_free(WFF0matrix);
  return 0;
}