#ifndef _UTIL_UTIL_H
#define _UTIL_UTIL_H

#include "constants.h"
#include "../NR/mdet.h"
#include "../NR/minverse.h"
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>

namespace Util {

/*************** handy ****************/
template<class T> void exch(T& a, T& b);

template <class T> int sign(const T& a);

template<class T> bool isnan(const T& x);

template<class T> bool isinf(const T& x);

template<class T> T min(const T& a, const T& b);

template<class T> T max(const T& a, const T& b);

template <class T> T abs(const T& a);

template<class T> T square(const T& a);

template<class T, class U> T cycle1D(const T& k, const U& N);

template<class T> 
bool readWord(const std::string& _w, T& _val);

template<class T, class U> 
bool readWord(const std::string& _w, T& _val, const U& _default);

template<class T>
bool readKeyItem(std::istream& _fin, const std::string& _keyword, T& _val);

template<class T, unsigned int N>
class Tuple {
public:
  T& operator()(unsigned int _k);
  const T& operator()(unsigned int _k) const;

private:
  T _data[N];
};

/**************** vector operations ****************/
template<class T> void vexch(T* a, T* b, int len);

template<class T> void vassign(const T& a, T* p, int len);

template<class T> bool vequal(const T& a, T* p, int len);

template<class T> void vcopy(const T* src, T* des, int len);

template<class T> void vcopy(const T& k, const T* src, T* des, int len);

template<class T> void vadd(const T* a, const T* b, T* c, int len);

template<class T> void vsub(const T* a, const T* b, T* c, int len);

template<class T> T vdotv(const T* a, const T* b, int len);

template<class T> T norm1(const T* a, int len);

template<class T> T norm2(const T* a, int len);

template<class T> void normalize(T* a, int len);

template<class T> void negv(T* a, int len);

template<class T> void kdotv(const T& k, const T* a, T* b, int len);

template<class T> T vvdist2(const T* a, const T* b, int len);

template<class T> void kdotvadd(const T& ka, const T* a, 
				const T& kb, const T* b, 
				T* c, int len);

template<class T> void kdotvadd3(const T& ka, const T* a, 
				 const T& kb, const T* b,
				 const T& kc, const T* c,
				 T* v, int len);

template<class T> void kdotvadd4(const T& ka, const T* a, 
				 const T& kb, const T* b,
				 const T& kc, const T* c,
				 const T& kd, const T* d,
				 T* v, int len);

template<class T>
void project(const T* a, const T* n, T* v, int len);

template<class T>
std::ostream& vdump(std::ostream& fout, const T* a, int len);

/**************** matrix operations ****************/
template<class T>
void mcopy(const T* mA, T* mX, int r, int s);

template<class T>
void mtrans(const T* mA, T* mX, int r, int s);

template<class T>
void madd(const T* mA, const T* mB, T* mX, int r, int s);

template<class T>
void msub(const T* mA, const T* mB, T* mX, int r, int s);

template<class T>
void kdotm(const T& k, const T* mA, T* mX, int r, int s);

template<class T>
void mdotv(const T* mA, const T* b, T* x, int r, int s);

template<class T>
void mdotm(const T* mA, const T* mB, T* mX, int r, int s, int t);

template<class T>
T mtrace(const T* mA, int n);

template<class T> 
T mdet(const T* mA, int n);

template<class T> 
void minverse(const T* mA, T* mX, int n);

/**************** statistics ****************/
template<class T> T sum(const T* v, int len);

template<class T> T mean(const T* v, int len);

template<class T> T var(const T* v, int len);

template<class T> void mean_var(const T* v, int len, T& _mean, T& _var);

/******************** series ********************/
template<class F> typename F::result_type kahan_sum_series(F& func, int bits);

template<class T> T expm1(T x);

/**************** misc ****************/
//a x b = c
template<class T> void cross3(const T* a, const T* b, T* c);

//(a, b, c)
template<class T> T cross_dot3(const T* a, const T* b, const T* c);

template<class T> void dual_axes(const T* a, const T* b, const T* c,
				 T* as, T* bs, T* cs);

template<class T> void build_axes(int n, T& ax);

template<class T> void linspace(const T& a, const T& z, T* result, int len);

template<class T> void cart2pol(const T& x, const T& y, T& r, T& th);

template<class T> void pol2cart(const T& r, const T& th, T& x, T& y);

template<class T> void rot(const T& th, T& x, T& y);

template<class T> void rot(const T* axis, const T& th, T* r);

/****************************************************************/
template<class T>
void exch(T& a, T& b) { T c(a); a = b; b = c; }

template <class T>
int sign(const T& a) { 
  if(a > 0) return 1; 
  else if(a < 0) return -1;
  else return 0;
}

template<class T>
bool isnan(const T& x) {
  return (x) != (x);
}

template<class T>
bool isinf(const T& x) {
  return 
    std::numeric_limits<T>::has_infinity &&
    x == std::numeric_limits<T>::infinity();
}

template<class T>
T min(const T& a, const T& b) { return ((a) < (b)) ? (a) : (b); }

template<class T>
T max(const T& a, const T& b) { return ((a) > (b)) ? (a) : (b); }

template <class T>
T abs(const T& a) { return ((a) < T(0) ? (-(a)) : (a)); }

template<class T>
T square(const T& a) { return a*a; }

template<class T, class U>
T cycle1D(const T& k, const U& N) { 
  return k+T(ceil(-k/double(N)))*N; 
}

template<class T> 
bool readWord(const std::string& _w, T& _val) {
  T t; std::istringstream iss;
  iss.clear();
  iss.str(_w);
  iss >> t;
  if(!iss.fail()) _val = t;
  return !iss.fail();
}

template<class T, class U> 
bool readWord(const std::string& _w, T& _val, const U& _default) {
  if(_w.empty()) { _val = _default; return false; }
  std::istringstream iss;
  iss.clear();
  iss.str(_w);
  iss >> _val;
  if(iss.fail()) _val = _default;
  return !iss.fail();
}

template<class T>
bool readKeyItem(std::istream& _fin, const std::string& _keyword, T& _val) {
  std::string _tword = "";
  _fin >> _tword;
  if(_tword == _keyword) _fin >> _val;
  else return false;
  return !_fin.fail();
}

template<class T, unsigned int N>
T& Tuple<T,N>::operator()(unsigned int _k) { return _data[_k]; } 

template<class T, unsigned int N>
const T& Tuple<T,N>::operator()(unsigned int _k) const { return _data[_k]; }

template<class T> 
void vexch(T* a, T* b, int len) {
  for(int j = 0; j < len; j++) exch(a[j],b[j]);
}

template<class T>
void vassign(const T& a, T* p, int len) { 
  for(int j = 0; j < len; j++) p[j] = a; }

template<class T> 
bool vequal(const T& a, T* p, int len) {
  for(int j = 0; j < len; j++) if(p[j] != a) return false; return true; }

template<class T>
void vcopy(const T* src, T* des, int len) { 
  for(int j = 0; j < len; j++) des[j] = src[j]; }

template<class T>
void vcopy(const T& k, const T* src, T* des, int len) { 
  for(int j = 0; j < len; j++) des[j] = k*src[j]; }

template<class T>
void vadd(const T* a, const T* b, T* c, int len) {
  for(int j = 0; j < len; j++) c[j] = a[j] + b[j];
}

template<class T>
void vsub(const T* a, const T* b, T* c, int len) {
  for(int j = 0; j < len; j++) c[j] = a[j]-b[j];
}

template<class T>
T vdotv(const T* a, const T* b, int len) {
  T v(0);
  for(int j = 0; j < len; j++) v += a[j]*b[j];
  return v;
}

template<class T>
T norm1(const T* a, int len) {
  T max_a(abs(a[0]));
  for(int j = 1; j < len; j++) 
    if(max_a < abs(a[j])) max_a = abs(a[j]);
  return max_a;
}

template<class T>
T norm2(const T* a, int len) {
  return T(std::sqrt(vdotv(a, a, len)));
}

template<class T>
void normalize(T* a, int len) {
  T norm_a(norm2(a,len));
  if(norm_a < EPS_DOUBLE) return;
  for(int j = 0; j < len; j++) a[j] /= norm_a;
}

template<class T>
void negv(T* a, int len) {
  for(int j = 0; j < len; j++) a[j] = -a[j];
}

template<class T>
void kdotv(const T& k, const T* a, T* b, int len) {
  for(int j = 0; j < len; j++) b[j] = a[j]*k;
}

template<class T>
T vvdist2(const T* a, const T* b, int len) {
  T v(0);
  for(int j = 0; j < len; j++) v += (a[j]-b[j])*(a[j]-b[j]);
  return v;
}

template<class T>
void kdotvadd(const T& ka, const T* a, const T& kb, const T* b, 
	      T* c, int len) {
  for(int j = 0; j < len; j++) c[j] = ka*a[j] + kb*b[j];
}

template<class T>
void kdotvadd3(const T& ka, const T* a, 
	       const T& kb, const T* b,
	       const T& kc, const T* c,
	       T* v, int len) {
  for(int j = 0; j < len; j++) v[j] = ka*a[j] + kb*b[j]+kc*c[j];
}

template<class T>
void kdotvadd4(const T& ka, const T* a, 
	       const T& kb, const T* b,
	       const T& kc, const T* c,
	       const T& kd, const T* d,
	       T* v, int len) {
  for(int j = 0; j < len; j++) v[j] = ka*a[j] + kb*b[j]+kc*c[j]+kd*d[j];
}

template<class T>
void project(const T* a, const T* n, T* v, int len) {
  T coef(vdotv(a,n,len)/vdotv(n,n,len));
  for(int j = 0; j < len; j++) v[j] = coef*n[j];
}

template<class T>
std::ostream& vdump(std::ostream& fout, const T* a, int len) {
  for(int j = 0; j < len; j++) fout << a[j] << (j==len-1?"":" ");
  return fout;
}

template<class T>
void mcopy(const T* mA, T* mX, int r, int s) {
  for(int j = 0, jmax = r*s; j < jmax; j++) mX[j] = mA[j];
}

template<class T>
void mtrans(const T* mA, T* mX, int r, int s) {
  int i, j;
  for(i = 0; i < r; i++)
    for(j = 0; j < s; j++) mX[j+i*s] = mA[i+j*r];
}

template<class T>
void madd(const T* mA, const T* mB, T* mX, int r, int s) {
  for(int j = 0, jmax = r*s; j < jmax; j++) mX[j] = mA[j]+mB[j];
}

template<class T>
void msub(const T* mA, const T* mB, T* mX, int r, int s) {
  for(int j = 0, jmax = r*s; j < jmax; j++) mX[j] = mA[j]-mB[j];
}

template<class T>
void kdotm(const T& k, const T* mA, T* mX, int r, int s) {
  for(int j = 0, jmax = r*s; j < jmax; j++) mX[j] = mA[j]*k;
}

template<class T>
void mdotv(const T* mA, const T* b, T* x, int r, int s) {
  int i, j;
  for(i = 0; i < r; i++) {
    x[i] = 0.0;
    for(j = 0; j < s; j++)
      x[i] += mA[i+j*r]*b[j];
  }
}

template<class T>
void mdotm(const T* mA, const T* mB, T* mX, int r, int s, int t) {
  T* entry;
  int i, j, k;
  for(i = 0; i < r; i++)
    for(j = 0; j < t; j++) {
      entry = &mX[i+j*r];
      for(k = 0, (*entry) = 0.0; k < s; k++) 
	(*entry) += mA[i+k*r]*mB[k+j*s];
    }
}

template<class T>
T mtrace(const T* mA, int n) {
  T t(0);
  for(int i = 0; i < n; i++) t += mA[i+i*n];
  return t;
}

template<class T> 
T mdet(const T* mA, int n) {
  int i, j;
  NR::Matrix<T> mAtmp(1,n,1,n);
  T** a = mAtmp.pointer();
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++) a[i+1][j+1] = mA[i+j*n];
  return NR::mdet(a,n);
}

template<class T> 
void minverse(const T* mA, T* mX, int n) {
  T d, *col, **a;
  int i, j, *indx;

  NR::Matrix<T> mAtmp(1,n,1,n);
  NR::Vector<T> vcol(1,n);
  NR::Vector<int> vindx(1,n);
  a = mAtmp.pointer();
  col=vcol.pointer();
  indx=vindx.pointer();
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++) a[i+1][j+1] = mA[i+j*n];
  NR::ludcmp(a,n,indx,&d);
  for(j = 1; j <= n; j++) {
    for(i = 1; i <= n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,n,indx,col);
    for(i = 1; i <= n; i++) mX[(i-1)+(j-1)*n] = col[i];
   } 
}

template<class T>
T sum(const T* v, int len) {
  T _sum(v[0]);
  for(int j = 1; j < len; j++) _sum += v[j];
  return _sum;
}

template<class T>
T mean(const T* v, int len) {
  return sum(v, len)/len;
}

template<class T>
T var(const T* v, int len) {
  T _mean, _var;
  mean_var(v,len,_mean,_var);
  return _var;
}

template<class T> 
void mean_var(const T* v, int len, T& _mean, T& _var) {
  /* provisional means algorithm */
  T tmp;
  _mean = v[0];
  _var = T(0);
  for(int k = 1; k < len; k++) {
    tmp = _mean;
    _mean = v[k]/(k+1)+(k*(_mean))/(k+1);
    _var = (v[k]-(_mean))*(v[k]-tmp)/k+((k-1)*(_var))/k;
  }
}

template<class F> 
typename F::result_type kahan_sum_series(F& func, int bits) {
  typedef typename F::result_type T;
  T factor = std::pow(T(2), bits);
  T result = func();
  T next_term, y, t;
  T carry = T(0);
  do {
    next_term = func();
    y = next_term-carry;
    t = result+y;
    carry = t-result;
    carry -= y;
    result = t;
  }while(std::fabs(result) < std::fabs(factor*next_term));
  return result;
}

//exp(x)-1 = sum_k(x^k/k!)
template<class T>
class expm1_series {
public:
  typedef T result_type;
  expm1_series(T x): k(0), m_x(x), m_term(1) {}
  T operator()() {
    k++;
    m_term *= m_x;
    m_term /= k;
    return m_term;
  }

private:
  int k;
  const T m_x;
  T m_term;
  expm1_series(const expm1_series&);
  expm1_series& operator= (const expm1_series&);
};

template<class T> 
T expm1(T x) {
  T a = T(std::fabs(x));
  if(a > T(0.5L))
    return std::exp(x)-T(1);
  if(a < std::numeric_limits<T>::epsilon())
    return x;
  expm1_series<T> s(x);
  return kahan_sum_series(s, std::numeric_limits<T>::digits+2);
}

template<class T>
void cross3(const T* a, const T* b, T* c) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = b[0]*a[2]-b[2]*a[0];
  c[2] = a[0]*b[1]-b[0]*a[1];
}

template<class T>
T cross_dot3(const T* a, const T* b, const T* c) {
  return
    (-a[2]*b[1]+a[1]*b[2])*c[0]+ 
    (+a[2]*b[0]-a[0]*b[2])*c[1]+
    (-a[1]*b[0]+a[0]*b[1])*c[2];
}

template<class T> void dual_axes(const T* a, const T* b, const T* c,
				 T* as, T* bs, T* cs) {
  T mat[9], invmat[9];
  vcopy(a, mat  , 3);
  vcopy(b, mat+3, 3);
  vcopy(c, mat+6, 3);
  minverse(mat, invmat, 3);
  for(int k = 0; k < 3; k++) as[k] = invmat[  3*k];
  for(int k = 0; k < 3; k++) bs[k] = invmat[1+3*k];
  for(int k = 0; k < 3; k++) cs[k] = invmat[2+3*k];
}

template<class T> 
void build_axes(int n, T& ax) {
  if(n == 3) return;
  if(n == 2) {
    cross3(ax[0], ax[1], ax[2]);
    normalize(ax[2],3);
  }
  if(n == 1) {
    ax[1][0] = ax[0][1]-ax[0][2]; 
    ax[1][1] = ax[0][2]-ax[0][0]; 
    ax[1][2] = ax[0][0]-ax[0][1];
    if(double(norm1(ax[1],3)) < EPS_DOUBLE) {
      ax[1][0] = ax[0][0]; 
      ax[1][1] = -ax[0][0]*ax[0][0]/ax[0][1]; 
      ax[1][2] = 0.0;
    }
    cross3(ax[0], ax[1], ax[2]);
    normalize(ax[1],3);
    normalize(ax[2],3);
  }
  if(n == 0) {
    ax[0][0] = 1.0; ax[0][1] = 0.0; ax[0][2] = 0.0;
    ax[1][0] = 0.0; ax[1][1] = 1.0; ax[1][2] = 0.0;
    ax[2][0] = 0.0; ax[2][1] = 0.0; ax[2][2] = 1.0;
  }
}

template<class T>
void linspace(const T& a, const T& z, T* result, int len) {
  for(int j = 0; j < len; j++) 
    result[j] = a+(z-a)*j/T(len-1);
}

template<class T>
void cart2pol(const T& x, const T& y, T& r, T& th) {
  r = std::sqrt(x*x+y*y);
  th = std::atan2(y, x);
}

template<class T>
void pol2cart(const T& r, const T& th, T& x, T& y) {
  x = r*std::cos(th);
  y = r*std::sin(th);
}

template<class T>
void rot(const T& th, T& x, T& y) {
  T tx = T(x*std::cos(th)-y*std::sin(th));
  T ty = T(x*std::sin(th)+y*std::cos(th));
  x = tx;
  y = ty;
}

template<class T>
void rot(const T* axis, const T& th, T* r) {
  T z[3], y[3], x[3], rz[3];
  vcopy(axis, z, 3);
  normalize(z, 3);
  cross3(z, r, y);
  T coef(vdotv(r,z,3));
  for(int j = 0; j < 3; j++) {
    rz[j] = coef*z[j];
    x[j] = r[j]-rz[j];
  }
  kdotvadd3<T>(std::cos(th),x,std::sin(th),y,1,rz,r,3);
}

}/* Util */

#endif
