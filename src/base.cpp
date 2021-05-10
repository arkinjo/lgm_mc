#include "base.hpp"
double log_add(const double a, const double b)
{
  /*
    log(c) = log(exp(a) + exp(b))
   */
  if(a > b) {
    return (a + log1p(exp(b - a)));
  } else {
    return (b + log1p(exp(a - b)));
  }
}

double log_dot(darray1& a, darray1& b, const int n)
{
  /*
    calculate 
    log(\sum_{k}[a_{k}*b_{k}])
    from [log(a_k)] and [log(b_k)].
   */
  darray1 s(boost::extents[n]); 
  double mv = neg_infinity;
  int imax = 0;
  for(int i = 0; i < n; ++i) {
    s[i] = a[i] + b[i];
    if(s[i] > mv) {
      mv = s[i];
      imax = i;
    }
  }
  double esum = 0.0;
  if(mv > neg_infinity && mv < pos_infinity) {
    for(int i = 0; i < n; ++i) {
      if(i != imax)
	esum += exp(s[i] - mv);
    }
  }

  double v = mv + log1p(esum);
  return v;
}

void log_mat_mult(darray2& a, darray2& b, darray2& c, const int l, const int m, const int n)
{
  /*
    calculate log(C) = log(AB) from log(A) and log(B).
    A: l X m
    B: m X n
    C: l X n
   */
  darray1 v(boost::extents[m]),w(boost::extents[m]);

  for(int i = 0; i < l; ++i) {
    for(int k = 0; k < m; ++k) {
      v[k] = a[i][k];
    }
    for(int j = 0; j < n; ++j) {
      for(int k = 0; k < m; ++k) {
	w[k] = b[k][j];
      }
      c[i][j] = log_dot(v,w,m);
    }
  }

  return;
}
