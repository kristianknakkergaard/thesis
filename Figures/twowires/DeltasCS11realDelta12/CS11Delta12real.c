#include <stdio.h>
#include <math.h>


int
main (void)
{
  double d_low = 0.71, d_up = 0.78, dd = 0.0001, dc = 0.7483;
  for (double d = d_low; d < d_up; d += dd)
  {
  if (d <= dc)
    {
    printf("%lg \t %lg\n", d, 0.0);
    }
  else
    printf("%lg \t %lg\n", d, 1.0);
  }
  

  return 0;
}