#include <stdio.h>
#include <math.h>


int
main (void)
{
  double d_low = 0.74, d_up = 0.775, dd = 0.000001, dc = 0.7586;
  for (double d = d_low; d < d_up; d += dd)
  {
  if (d < dc)
    {
    printf("%lg \t %lg\n", d, 0.0);
    }
  else
    printf("%lg \t %lg\n", d, 1.0);
  }
  

  return 0;
}