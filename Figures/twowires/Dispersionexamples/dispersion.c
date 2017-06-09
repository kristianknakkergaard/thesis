#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>


double EFplus(double epsilon, double D11, double D12)
{
	return sqrt(epsilon * epsilon + (D11 + D12) * (D11 + D12));
}

double EFminus(double epsilon, double D11, double D12)
{
  return sqrt(epsilon * epsilon + (D11 - D12) * (D11 - D12));
}
double EF(double epsilon, double D11, double D12)
{
  return sqrt(epsilon * epsilon + D11 * D11 + D12 * D12);
}


int
main (void)
{
  double k_low = -5.0, k_up = 5.0, dk = 0.001;
 
  double epsilonk;
  double D11k;
  double D12k;

  printf("%s \t %s \t %s \t %s \t %s \n", "k", "Efreegas", "EFkplus", "EFkminus", "EFk");

  for (double k = k_low; k < k_up; k += dk)
  {
    epsilonk 	= k * k - 1.0;
    D11k 		= 0.3 * k;
    D12k 		= 0.1;
    printf("%lg \t %lg \t %lg \t %lg \t %lg \n", k, fabs(epsilonk), EF(epsilonk, D11k, D12k), EFplus(epsilonk, D11k, D12k), EFminus(epsilonk, D11k, D12k) );
  }
  
  return 0;
}