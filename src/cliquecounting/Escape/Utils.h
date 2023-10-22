#ifndef ESCAPE_UTILS_H_
#define ESCAPE_UTILS_H_

#include <algorithm>

namespace Escape
{

//x choose n.
template <int n, typename T>
T choose(T x)
{
  T num = 1;
  T den = x;
  for (int i = 1; i < n; ++i)
  {
    num *= (i + 1); //this can easily overflow!!!
    den *= (x - i);
  }
  
  return den / num;
}


//max(0, x - y)
//Safe to use with unsigned types
template <typename T>
T czsub(T x, T y)
{
  return x > y ? x - y : 0;
}



}
#endif
