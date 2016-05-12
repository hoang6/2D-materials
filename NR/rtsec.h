#ifndef _NR_RTSEC_H
#define _NR_RTSEC_H

#include "util.h"

namespace NR {

/************************************************************/
template<class F, class T>
void rtsec(T x1, T x2, T tol, F& func, int& iter, T& root);

class Error_rtsec {};

/************************************************************/
template<class F, class T>
void rtsec(T x1, T x2, T tol, F& func, int& iter, T& root) {
  T fl,f,dx,swap,xl;

  fl = func.func(x1);
  f = func.func(x2);
  if(ABS(fl) < ABS(f)) {
    root = x1;
    xl = x2;
    swap = fl;
    fl = f;
    f = swap;
  } else {
    xl = x1;
    root = x2;
  }
  for(iter = 1;iter <= 100; iter++) {
    dx = (xl-root)*f/(f-fl);
    xl = root;
    fl = f;
    root += dx;
    f = func.func(root);
    if(ABS(dx) < tol || f == 0.0) return;
  }
  throw Error_rtsec();
}

} /* NR */

#endif
