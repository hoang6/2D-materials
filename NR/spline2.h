#ifndef _NR_SPLINE2_H
#define _NR_SPLINE2_H

#include "spline1.h"

namespace NR {

/********** definitions **********/
template<class T>
class Spline2 {
public:
  Spline2();
  Spline2(int _m, int _n);
  Spline2(T _xrow[], T _xcol[], T** _y, int _m, int _n);
  ~Spline2();
  Spline2(const Spline2& spline2);
  Spline2& operator= (const Spline2& spline2);
  void resize(int _m, int _n);
  void train(T _xrow[], T _xcol[], T** _y);
  void train(T __xrow, T __xcol, T __y, int __m, int __n);
  void prepare();
  T evaluate(T __xrow, T __xcol);
  
protected:
  void reset();
  void allocate(int _m, int _n);
  void deallocate();
  void copy(T _xrow[], T _xcol[], T** _y, int m0, int n0);

  T *xrow, *xcol, **y, **y2;
  int m, n;
  bool allocated;
};

template<class T>
void splie2(T x1a[], T x2a[], T **ya, int m, int n, T **y2a);

template<class T>
void splin2(T x1a[], T x2a[], T **ya, T **y2a, int m, int n,
	    T x1, T x2, T *y);

/********** definitions **********/
template<class T>
Spline2<T>::Spline2() { reset(); }

template<class T>
Spline2<T>::Spline2(int _m, int _n) { reset(); allocate(_m, _n); }

template<class T>
Spline2<T>::Spline2(T _xrow[], T _xcol[], T** _y, int _m, int _n) {
  reset();
  allocate(_m,_n);
  copy(_xrow, _xcol, _y, 0, 0);
}

template<class T>
Spline2<T>::~Spline2() {
  if(allocated) deallocate();
}

template<class T>
Spline2<T>::Spline2(const Spline2& spline2) {
  reset();
  *this = spline2;
}

template<class T>
Spline2<T>& Spline2<T>::operator= (const Spline2<T>& spline2) {
  if(this != &spline2) {
    if(!spline2.allocated && !allocated) return;
    else if(!spline2.allocated && allocated) { deallocate(); }
    else {
      resize(spline2.m, spline2.n);
      copy(spline2.xrow, spline2.xcol, spline2.y, 1, 1);      
    }
  }
  return *this;
}

template<class T>
void Spline2<T>::resize(int _m, int _n) {
  if(m == _m && n == _n) return;
  if(allocated) deallocate();
  allocate(_m,_n);
}

template<class T>
void Spline2<T>::train(T _xrow[], T _xcol[], T** _y) {
  copy(_xrow,_xcol,_y,0,0);
}

template<class T>
void Spline2<T>::train(T __xrow, T __xcol, T __y, int __m, int __n) {
  xrow[__m+1] = __xrow;
  xcol[__n+1] = __xcol;
  y[__m+1][__n+1] = __y;
}

template<class T>
inline void Spline2<T>::prepare() {
  splie2(xrow,xcol,y,m,n,y2);
}

template<class T>
inline T Spline2<T>::evaluate(T __xrow, T __xcol) {
  T __y;
  splin2(xrow, xcol, y, y2, m, n, __xrow, __xcol, &__y);
  return __y;
}

template<class T>
void Spline2<T>::reset() {
  xrow = xcol = NULL;
  y = y2 = NULL;
  m = n = 0;
  allocated = false;
}

template<class T>
void Spline2<T>::allocate(int _m, int _n) {
  xrow = vector<T>(1,_m);
  xcol = vector<T>(1,_n);
  y = matrix<T>(1,_m,1,_n);
  y2 = matrix<T>(1,_m,1,_n);
  m = _m; n = _n;
  allocated = true;
}

template<class T>
void Spline2<T>::deallocate() {
  free_vector<T>(xrow,1,m);
  free_vector<T>(xcol,1,m);
  free_matrix<T>(y,1,m,1,n);
  free_matrix<T>(y2,1,m,1,n);
  reset();
}

template<class T>
void Spline2<T>::copy(T _xrow[], T _xcol[], T** _y, int m0, int n0) {
  int r, c;
  for(r = 1; r <= m; r++) xrow[r] = _xrow[r-1+m0];
  for(c = 1; c <= n; c++) xcol[c] = _xcol[c-1+n0];
  for(r = 1; r <= m; r++)
    for(c = 1; c <= n; c++) y[r][c] = _y[r-1+m0][c-1+n0];
}

template<class T>
void splie2(T x1a[], T x2a[], T **ya, int m, int n, T **y2a) {
  int j;
  for (j=1;j<=m;j++)
    spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
}

template<class T>
void splin2(T x1a[], T x2a[], T **ya, T **y2a, int m, int n,
	    T x1, T x2, T *y) {
  int j;
  T *ytmp,*yytmp;
  
  Vector<T> vytmp(1,m);
  Vector<T> vyytmp(1,m);
  ytmp=vytmp.pointer();
  yytmp=vyytmp.pointer();
  for (j=1;j<=m;j++)
    splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  splint(x1a,yytmp,ytmp,m,x1,y);
}

} /* NR */

#endif
