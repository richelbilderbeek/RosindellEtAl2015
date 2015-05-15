#include "helper.h"

#include <cmath>

std::string alpha_format(long x)
{
  const char c = 'A' + x;
  return std::string(1,c);
}


std::string alpha_format(long x_in, int num_digits)
{
  std::string alpha_res;
  alpha_res = "";
  long x;
  x = x_in;
  if (pow(26.0,double(num_digits)) >= x)
  {
    for (int i = 0 ; i < num_digits ; ++i)
    {
      long remainder = (x % 26);
      alpha_res = alpha_format(remainder) + alpha_res;
      x -= remainder;
      x = x/26;
    }
  }
  else
  {
    alpha_res = "_ERROR_";
  }
  return alpha_res;
}

std::string alpha_format_auto(long x_in)
{
  long upto = 25;
  long num_digits = 1;
  long x = x_in;
  while (upto < x_in)
  {
    x -= long(floor(pow(26.0,num_digits)));
    num_digits ++;
    upto += long(floor(pow(26.0,num_digits)));
  }
  return (alpha_format(x,num_digits));
}
