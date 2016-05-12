#ifndef _UTIL_CONJGRAD_H
#define _UTIL_CONJGRAD_H

#include "MinimizerBasic.h"
#include "Accessors.h"
#include "util.h"

namespace Util {

template<class F, class T = double> class ConjGrad;

template<class F, class T>
bool conjGrad(T x[], int n, T tol, int maxIterNum, F& funcd, 
	      int& iterNum, T& fmin);

/************************************************************/
template<class F, class T>
class ConjGrad: public MinimizerBasic<F, T> {
public:
  ConjGrad();
  ConjGrad(T _tol, int _maxIterNum);
  ConjGrad(const ConjGrad& conjGrad);
  virtual ~ConjGrad();
  ConjGrad& operator= (const ConjGrad& conjGrad);
  virtual void reset();
  virtual bool minimize(T* _x, int _n, F& _funcd);
  virtual typename MinimizerBasic<F, T>::Style style() const;
  T minFunctionValue() const;

protected:
  void initialize(T* _x, int _n, F& _funcd);
  void finalize(int _iter);
  void shift(T& a, T& b, T& c, const T& d);
  void bracket_min(T* ax, T* bx, T* cx, T* fa, T* fb, T* fc);
  T absign(const T& a, const T& b);
  T projection(T x);
  T parab_intrp(T ax, T bx, T cx, T tol, T* xmin);
  void minimize_1d(T& fret);

  T* p;
  int n;
  F* pfuncd;
  T* g;
  T* h;
  T* xi;
  T* xt;
  T fmin;
  T max_dx;
  int error_flag;

private:
  void _reset();
};

/************************************************************/
template<class F, class T>
inline ConjGrad<F, T>::ConjGrad() { _reset(); }

template<class F, class T>
inline ConjGrad<F, T>::ConjGrad(T _tol, int _maxIterNum): MinimizerBasic<F, T>(_tol, _maxIterNum, 0.0) {}

template<class F, class T>
inline ConjGrad<F, T>::ConjGrad(const ConjGrad<F, T>& conjGrad) { *this = conjGrad; }

template<class F, class T>
inline ConjGrad<F, T>::~ConjGrad() {}

template<class F, class T>
ConjGrad<F, T>& ConjGrad<F, T>::operator= (const ConjGrad<F, T>& conjGrad) {
  if(this != &conjGrad) {
    this->MinimizerBasic<F, T>::operator= (conjGrad);
  }
  return *this;
}

template<class F, class T>
inline void ConjGrad<F, T>::reset() { _reset(); }

template<class F, class T>
bool ConjGrad<F, T>::minimize(T* _x, int _n, F& _funcd) {
  int j, iter;
  T gg, gam, fp, dgg;
  T max_p, max_g;

  initialize(_x, _n, _funcd);
  iter = 0;
  error_flag = 0;
  fp = (*pfuncd)(p);
  (*pfuncd)(p, xi);

  max_p = 0.0;
  max_g = 0.0;
  for(j = 0; j < n; j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
    max_p = max(max_p,abs(p[j]));
    max_g = max(max_g,abs(g[j]));
  }
  if(max_g == 0.0) {
    finalize(iter);
    return true;
  }
  if((max_dx = 100.0*max_p/max_g) < 100.0) max_dx = 100.0;
  //std::cout << "ConjGrad::minimize max_dx = " << max_dx << std::endl;

  for(iter = 1; iter <= this->maxIterNumber(); iter++) {
    //debug
    (*pfuncd)(p, g, fp);

    if(this->isbad(fp)) {
      finalize(iter);
      return false;
    }
    
    minimize_1d(fmin);
    if(error_flag) {
      for(j = 0; j < n; j++) {
	xi[j]=h[j]=g[j];
      }
      error_flag = 0;
#ifndef NR_SILENT
      std::cerr << "Search directions reset in ConjGrad::minimize" << std::endl;
#endif
      continue;
    }

    if(2.0*abs(fmin-fp) <= this->tolerance()*(abs(fmin)+abs(fp)+1.0e-10)) {
      finalize(iter);
      return true;
    }
    
    fp = (*pfuncd)(p);
    (*pfuncd)(p, xi);
    dgg = gg = 0.0;
    for(j = 0; j < n; j++) {
      gg += g[j]*g[j];
      dgg += (xi[j]+g[j])*xi[j];
    }
    if(gg == 0.0) {
      finalize(iter);
      return true;
    }
    gam = dgg/gg;
    for(j = 0; j < n; j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }//end for(iter = 1...

#ifndef NR_SILENT
  std::cerr << "Too many iterations in ConjGrad::minimize." << std::endl;
#endif

  finalize(iter-1);
  return false;
}

template<class F, class T>
inline typename MinimizerBasic<F, T>::Style ConjGrad<F, T>::style() const { return MinimizerBasic<F, T>::FUNC_DFUNC; }

template<class F, class T>
inline T ConjGrad<F, T>::minFunctionValue() const { return fmin; }

template<class F, class T>
void ConjGrad<F, T>::initialize(T* _x, int _n, F& _funcd) {
  p = _x;
  n = _n;
  pfuncd = &_funcd;
  g = new T[n];
  h = new T[n];
  xi = new T[n];
  xt = new T[n];
}

template<class F, class T>
void ConjGrad<F, T>::finalize(int _iter) {
  delete [] xt;
  delete [] xi;
  delete [] h;
  delete [] g;

  this->iterNum = _iter;
}

template<class F, class T>
inline void ConjGrad<F, T>::shift(T& a, T& b, T& c, const T& d) { 
  a = b; b = c; c = d; 
}

template<class F, class T>
void ConjGrad<F, T>::bracket_min(T* ax, T* bx, T* cx, T* fa, T* fb, T* fc) {
  T ulim,u,r,q,fu,dum;
  *fa = projection(*ax);
  *fb = projection(*bx);
  if(*fb > *fa) {
    shift(dum,*ax,*bx,dum);
    shift(dum,*fb,*fa,dum);
  }
  *cx = (*bx)+1.618034*(*bx-*ax);
  *fc = projection(*cx);
  if(error_flag) return;
  while(*fb > *fc) {
    r = (*bx-*ax)*(*fb-*fc);
    q = (*bx-*cx)*(*fb-*fa);
    u = (*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*absign(max(abs(q-r),1.0e-20),q-r));
    ulim = (*bx)+100.0*(*cx-*bx);
    if((*bx-u)*(u-*cx) > 0.0) {
      fu = projection(u);
      if(error_flag) return;
      if (fu < *fc) {
	*ax = (*bx);
	*bx = u;
	*fa = (*fb);
	*fb = fu;
	return;
      }
      else if (fu > *fb) {
	*cx = u;
	*fc = fu;
	return;
      }
      u = (*cx)+1.618034*(*cx-*bx);
      fu = projection(u);
      if(error_flag) return;
    }
    else if ((*cx-u)*(u-ulim) > 0.0) {
      fu = projection(u);
      if(error_flag) return;
      if (fu < *fc) {
	shift(*bx,*cx,u,*cx+1.618034*(*cx-*bx));
	shift(*fb,*fc,fu, projection(u));
	if(error_flag) return;
      }
    }
    else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u = ulim;
      fu = projection(u);
      if(error_flag) return;
    }
    else {
      u = (*cx)+1.618034*(*cx-*bx);
      fu = projection(u);
      if(error_flag) return;
    }
    shift(*ax,*bx,*cx,u);
    shift(*fa,*fb,*fc,fu);
  }//end while(*fb > *fc)
}

template<class F, class T>
T ConjGrad<F, T>::absign(const T& a, const T& b) {
  return ((b) >= T(0) ? abs(a) : (-abs(a)));
}

template<class F, class T>
T ConjGrad<F, T>::projection(T x) {
  if(abs(x) >= max_dx) { 
    error_flag = 1; 
    return 0.0; 
  }
  for (int j = 0; j < n; j++) 
    xt[j] = p[j]+x*xi[j];
  return (*pfuncd)(xt);
}

template<class F, class T>
T ConjGrad<F, T>::parab_intrp(T ax, T bx, T cx, T tol, T* xmin) {
  int iter;
  T a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  T e = 0.0;
  d = 0.0;
  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = projection(x);
  if(error_flag) return 0.0;
  for(iter = 1; iter <= 100; iter++) {
    xm = 0.5*(a+b);
    tol2 = 2.0*(tol1 = tol*abs(x)+1.0e-10);
    if(abs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if(abs(e) > tol1) {
      r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if(q > 0.0) p = -p;
      q = abs(q);
      etemp = e;
      e = d;
      if(abs(p) >= abs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d = 0.3819660*(e=(x >= xm ? a-x : b-x));
      else {
	d = p/q;
	u = x+d;
	if (u-a < tol2 || b-u < tol2) d = absign(tol1,xm-x);
      }
    }
    else d = 0.3819660*(e=(x >= xm ? a-x : b-x));
    u = (abs(d) >= tol1 ? x+d : x+absign(tol1,d));
    fu = projection(u);
    if(error_flag) return 0.0;
    if(fu <= fx) {
      if(u >= x) a = x; 
      else b = x;
      shift(v,w,x,u);
      shift(fv,fw,fx,fu);
    }
    else {
      if (u < x) a = u; 
      else b = u;
      if (fu <= fw || w == x) {
	v = w;
	w = u;
	fv = fw;
	fw = fu;
      }
      else if (fu <= fv || v == x || v == w) {
	v = u;
	fv = fu;
      }
    }
  }//end for(iter = 1...
#ifndef NR_SILENT
  std::cerr << "Too many iterations in ConjGrad::parab_intrp." << std::endl;
#endif
  *xmin=x;
  return fx;
}

template<class F, class T>
void ConjGrad<F, T>::minimize_1d(T& fret) {
  int j;
  T xx,xmin,fx,fb,fa,bx,ax;
  ax = 0.0;
  xx = 1.0;
  bracket_min(&ax,&xx,&bx,&fa,&fx,&fb);
  if(error_flag) return;
  xmin = T(0);
  fret = parab_intrp(ax,xx,bx,2.0e-4,&xmin);
  if(error_flag) return;
  for(j = 0; j < n; j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
}

template<class F, class T>
void ConjGrad<F, T>::_reset() {
  this->tolerance(1.0e-9);
  this->maxIterNumber(10000);
  this->deltaTime(0.0);
}

/************************************************************/
template<class F, class T>
bool conjGrad(T x[], int n, T tol, int maxIterNum, F& funcd, 
	      int& iterNum, T& fmin) {
  ConjGrad<F, T> cg(tol, maxIterNum);
  bool flag = cg.minimize(x, n, funcd);
  iterNum = cg.iterNumber();
  fmin = cg.minFunctionValue();
  return flag;
}

} /* Util */

#endif
