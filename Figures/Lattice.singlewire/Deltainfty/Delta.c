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

double Vindx (double x, double xi, double strength)
{
  return - strength * exp( - sqrt(2.0) * fabs(x) / xi );
}

double WFF0( double k, double kprime, double xi, double strength, int N )
{
  double W = 0.0;

  for (int i = 1; i < N; ++i)
  {
    double x  = (double) i;
    W += 1.0 / 2.0 * ( cos((kprime - k) * i ) - cos( (kprime + k) * i ) ) * Vindx( x, xi, strength ); 
  }

  return W; 
}

double Deltaguess(double k, double a, double b)
{
  return a * sin(k) + b * sin(2.0 * k);
}


int
main (void)
{
 /*restricted variables */
  double t1       = 0.0;  
  double strength = 4.0; /*total strength of interaction*/

  /*free parameters:*/
  double xi_low   = 1.0,  dxi   = 0.1; /*we will eventually want to control the range rather than the Bose gas parameter*/
  double fill_low = 0.0, dfill  = 0.01;
  double t2       = 1.0;
  int    N        = 100; 
  double Ndouble  = (double)N;
  int    Nxi      = 61;
  int    Nfill    = 101;
  
  /*k-values:*/
  double k_low = - M_PI, dk = 2.0 * M_PI / Ndouble;

  /*iterations and convergence: */

  int    iterations     = 200;
  double convergence    = 1e-4;
  double filltolerance  = 5e-4;

  /*function variables:*/
  double epsilonkprime;
  double EFkprime;
  double Deltakprime;
  
  double gapintegrand;
  double gapintegral;
  double D_maxk;

  double k;
  double kprime;
  int check1k = 0, check2k = 0, flag1k = 0, flag2k = 0;

  /*vectors and matrices:*/
  gsl_vector *D1k           = gsl_vector_calloc(N);
  gsl_vector *D2k           = gsl_vector_calloc(N);
  gsl_matrix *WFF0matrix    = gsl_matrix_calloc(N, N);

  /*for the chemical potential*/
  double mu_low     = -3.0, dmu = 0.001, mu1k, mu2k;
  int Nmu           = 6000;
  double muintegrand, muintegral = 0.0, fillingnow1k = 0.0, fillingnow2k = 0.0, fillingprior1k, fillingprior2k;
  double mu_min_value;
  
  double mu_guess;
  gsl_vector *mu_vector     = gsl_vector_calloc(Nmu);
  gsl_vector *find_min_mu   = gsl_vector_calloc(Nmu);
  gsl_vector *fillingvector = gsl_vector_calloc(Nmu);

  for (int imu = 0; imu < Nmu; ++imu)
  {
    mu_guess = (double) mu_low + imu * dmu;
    gsl_vector_set(mu_vector, imu, mu_guess);
  }

  /*variables for ground state energy*/
  double E0integrand, E01k, E02k;

  /*variables for winding */
  double winding1k, winding2k;
  double kprimenext;
  double epsilonkprimenext; 
  int Nhalf = (int) Ndouble / 2.0; 
  
  /* start of xi-fill double loop*/
  double xi, fill;

  printf("%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n", "fill", "xi", "winding1k", "winding2k", "E01k + fillingnow1k * mu", "E02k + fillingnow2k * mu", "Delta1kmax", "Delta2kmax" );

  for (int jxi = 0; jxi < Nxi; ++jxi)
  {
    xi = xi_low + (double)jxi * dxi;

    /*interaction:*/
    for (int i = 0; i < N; ++i)
    {
      k = ((double)i) * dk + k_low;

      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime = ((double)iprime) * dk + k_low;
        gsl_matrix_set( WFF0matrix, i, iprime, WFF0(k, kprime, xi, strength, N ) );
      }
    }

    for (int jfill = 0; jfill < Nfill; ++jfill)
    {
      fill = fill_low + (double)jfill * dfill;
      mu1k = -3.0 + 2.0 * 3.0 * fill, check1k = 0, check2k = 0, flag1k = 0, flag2k = 0;
      mu2k = -3.0 + 2.0 * 3.0 * fill;
      
      /*sin(1k)-like: */

      /*start guess for Delta:*/
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gsl_vector_set( D1k, i, Deltaguess(k, 1.0, 0.0) );
      }

      for (int j = 0; j < iterations; ++j)
      {
        D_maxk = gsl_vector_max(D1k);
        for (int i = 0; i < N; ++i)
        {
          k = ((double)i) * dk + k_low;
          gapintegral = 0.0;
          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime        = ((double)iprime) * dk + k_low;
            epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu1k;
            Deltakprime   = gsl_vector_get(D1k, iprime);
            EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
            gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
            gapintegral  += gapintegrand; 
          }
          gsl_vector_set(D1k, i, gapintegral);
        }
        for (int imu = 0; imu < Nmu; ++imu)
        {
          mu_guess = mu_low + (double)imu * dmu;
          muintegral = 0.0;
          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime        = ((double)iprime) * dk + k_low;
            epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu_guess;
            Deltakprime   = gsl_vector_get(D1k, iprime);
            EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
            muintegrand   = 1.0 / (2.0 * Ndouble) * (1.0 -  epsilonkprime / EFkprime);
            muintegral   += muintegrand;
          }
          gsl_vector_set(find_min_mu, imu, pow(fill - muintegral, 2.0) );
          gsl_vector_set(fillingvector, imu, muintegral );
        }
        mu_min_value  = gsl_vector_get(mu_vector,     gsl_vector_min_index(find_min_mu));
        fillingnow1k  = gsl_vector_get(fillingvector, gsl_vector_min_index(find_min_mu));

        if ( fabs(mu1k - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D1k) - D_maxk)/D_maxk < convergence && fabs( fillingprior1k - fillingnow1k ) < filltolerance && fabs( fill - fillingnow1k ) < 10 * filltolerance )
        {
          break;
        }

        if ( fabs(mu1k - mu_min_value) < convergence && fabs(mu1k) < convergence && fabs(gsl_vector_max(D1k) - D_maxk)/D_maxk < convergence && fabs( fillingprior1k - fillingnow1k ) < filltolerance && fabs( fill - fillingnow1k ) < 10 * filltolerance )
        {
          break;
        }
        if ( fabs(mu1k - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D1k)) < convergence && fabs( fillingprior1k - fillingnow1k ) < filltolerance && fabs( fill - fillingnow1k ) < 10 * filltolerance )
        {
          break;
        }

        if ( fabs(mu1k - mu_min_value) < convergence && fabs(mu1k) < convergence && fabs(gsl_vector_max(D1k)) < convergence && fabs( fillingprior1k - fillingnow1k ) < filltolerance && fabs( fill - fillingnow1k ) < 10 * filltolerance )
        {
          break;
        }

        fillingprior1k = fillingnow1k;
        mu1k           = mu_min_value;
        ++check1k; 
      }

      if (check1k == iterations)
      {
        flag1k = 1;
      }

      /*Ground state energy: */
      E01k = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu1k;
        Deltakprime   = gsl_vector_get(D1k, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

        E0integrand   = - 1.0 / (4.0 * Ndouble) * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
        E01k         += E0integrand; 
      }

      /*Topological invariant*/
      winding1k = 0.0; 
      for (int iprime = 0; iprime < Nhalf; ++iprime)
      {
        kprime            = ((double)iprime      ) * dk + k_low;
        kprimenext        = ((double)iprime + 1.0) * dk + k_low;

        epsilonkprime     = - t1 * cos(kprime)     - t2 * cos(2.0 * kprime)     - mu1k;
        epsilonkprimenext = - t1 * cos(kprimenext) - t2 * cos(2.0 * kprimenext) - mu1k;

        if (epsilonkprime < 0 && epsilonkprimenext > 0)
        {
          winding1k += copysign(1.0, gsl_vector_get(D1k, iprime)); 
        }
        if (epsilonkprime > 0 && epsilonkprimenext < 0)
        {
          winding1k -= copysign(1.0, gsl_vector_get(D1k, iprime)); 
        }
      }
      

      /*sin(2k)-like: */

      /*start guess for Delta:*/
      for (int i = 0; i < N; ++i)
      {
        k = ((double)i) * dk + k_low;
        gsl_vector_set( D2k, i, Deltaguess(k, 0.0, 1.0) );
      }

      for (int j = 0; j < iterations; ++j)
      {
        D_maxk = gsl_vector_max(D2k);
        for (int i = 0; i < N; ++i)
        {
          k = ((double)i) * dk + k_low;
          gapintegral = 0.0;
          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime        = ((double)iprime) * dk + k_low;
            epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu2k;
            Deltakprime   = gsl_vector_get(D2k, iprime);
            EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
            gapintegrand  = - 1.0 / Ndouble * gsl_matrix_get(WFF0matrix, i, iprime) * Deltakprime / (2.0 * EFkprime);
            gapintegral  += gapintegrand; 
          }
          gsl_vector_set(D2k, i, gapintegral);
        }
        for (int imu = 0; imu < Nmu; ++imu)
        {
          mu_guess = mu_low + (double)imu * dmu;
          muintegral = 0.0;
          for (int iprime = 0; iprime < N; ++iprime)
          {
            kprime        = ((double)iprime) * dk + k_low;
            epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu_guess;
            Deltakprime   = gsl_vector_get(D2k, iprime);
            EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);
            muintegrand   = 1.0 / (2.0 * Ndouble) * (1.0 -  epsilonkprime / EFkprime);
            muintegral   += muintegrand;
          }
          gsl_vector_set(find_min_mu, imu, pow(fill - muintegral, 2.0) );
          gsl_vector_set(fillingvector, imu, muintegral );
        }
        mu_min_value  = gsl_vector_get(mu_vector,     gsl_vector_min_index(find_min_mu));
        fillingnow2k  = gsl_vector_get(fillingvector, gsl_vector_min_index(find_min_mu));

        if ( fabs(mu2k - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D2k) - D_maxk)/D_maxk < convergence && fabs( fillingprior2k - fillingnow2k ) < filltolerance && fabs( fill - fillingnow2k ) < 10 * filltolerance )
        {
          break;
        }

        if ( fabs(mu2k - mu_min_value) < convergence && fabs(mu2k) < convergence && fabs(gsl_vector_max(D2k) - D_maxk)/D_maxk < convergence && fabs( fillingprior2k - fillingnow2k ) < filltolerance && fabs( fill - fillingnow2k ) < 10 * filltolerance )
        {
          break;
        }
        if ( fabs(mu2k - mu_min_value)/mu_min_value < convergence && fabs(gsl_vector_max(D2k)) < convergence && fabs( fillingprior2k - fillingnow2k ) < filltolerance && fabs( fill - fillingnow2k ) < 10 * filltolerance )
        {
          break;
        }

        if ( fabs(mu2k - mu_min_value) < convergence && fabs(mu2k) < convergence && fabs(gsl_vector_max(D2k)) < convergence && fabs( fillingprior2k - fillingnow2k ) < filltolerance && fabs( fill - fillingnow2k ) < 10 * filltolerance )
        {
          break;
        }
        
        fillingprior2k = fillingnow2k;
        mu2k           = mu_min_value;
        ++check2k; 
      }

      if (check2k == iterations )
      {
        flag2k = 1;
      }
      /*Ground state energy: */
      E02k = 0.0;
      for (int iprime = 0; iprime < N; ++iprime)
      {
        kprime        = ((double)iprime) * dk + k_low;
        epsilonkprime = - t1 * cos(kprime) - t2 * cos(2.0 * kprime) - mu2k;
        Deltakprime   = gsl_vector_get(D2k, iprime);
        EFkprime      = pow(epsilonkprime * epsilonkprime + Deltakprime * Deltakprime, 1.0 / 2.0);

        E0integrand   = - 1.0 / (4.0 * Ndouble) * pow(epsilonkprime - EFkprime, 2.0) / EFkprime;  
        E02k         += E0integrand; 
      }

      /*Topological invariant*/
      winding2k = 0.0; 
      for (int iprime = 0; iprime < Nhalf; ++iprime)
      {
        kprime            = ((double)iprime      ) * dk + k_low;
        kprimenext        = ((double)iprime + 1.0) * dk + k_low;

        epsilonkprime     = - t1 * cos(kprime)     - t2 * cos(2.0 * kprime)     - mu2k;
        epsilonkprimenext = - t1 * cos(kprimenext) - t2 * cos(2.0 * kprimenext) - mu2k;

        if (epsilonkprime < 0 && epsilonkprimenext > 0)
        {
          winding2k += copysign(1.0, gsl_vector_get(D2k, iprime)); 
        }
        if (epsilonkprime > 0 && epsilonkprimenext < 0)
        {
          winding2k -= copysign(1.0, gsl_vector_get(D2k, iprime)); 
        }
      }

      printf("%lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %i \t %i  \n", fill, xi, fabs(winding1k), fabs(winding2k), E01k + fillingnow1k * mu1k, E02k + fillingnow2k * mu2k, gsl_vector_max(D1k), gsl_vector_max(D2k), flag1k, flag2k);

      fprintf(stderr, "fill = %lg, \t xi = %lg, \t winding1k = %lg, \t winding2k = %lg, \t flag1k = %i, \t flag2k = %i \n", fill, xi, fabs(winding1k), fabs(winding2k), flag1k, flag2k );
    }
    printf("\n");
  }

 
  gsl_vector_free(mu_vector);
  gsl_vector_free(find_min_mu);
  gsl_vector_free(fillingvector);
  gsl_vector_free(D1k);
  gsl_vector_free(D2k);
  gsl_matrix_free(WFF0matrix);
  return 0;
}