#include <stdio.h>
#include <math.h>

double Delta1k(double k, double a)
{
  return a * sin(k);
}

double Delta2k(double k, double a)
{
  return a * sin(2.0 * k);
}
double epsilonk(double k, double t1, double t2, double mu)
{
  return - t1 * cos(k) - t2*cos(2.0 * k) - mu;
}


int
main (void)
{
  double t1 = 1.0; 
  double t2 = 1.0; 
  double dk = 0.00001;
  double mu1 = -0.5, mu2 = 0.5;

  printf("%s \t %s \t %s \t %s \t %s \n", "k", "epsilonk1", "epsilonk2", "Delta1k", "Delta2k" );
  for (double k = -M_PI; k < M_PI; k+=dk)
  {
    printf("%lg \t %lg \t %lg \t %lg \t %lg \n", k, epsilonk(k, t1, t2, mu1), epsilonk(k, t1, t2, mu2), Delta1k(k, 0.5), Delta2k(k, 0.5) );
  }
  return 0;
}