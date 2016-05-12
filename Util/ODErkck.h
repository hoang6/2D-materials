#ifndef _UTIL_ODERKCK_H
#define _UTIL_ODERKCK_H

#include "constants.h"
#include "Array1D.h"
#include "util.h"
#include <vector>
#include <cmath>

namespace Util {

/* 
 * Cash-Karp method: 
 * a 4/5 embedded R-K scheme
 * u'=f(t,u) ==> 
 *   F.d(t,u,du) ~ F.d(T, const Data&, Data&) 
 *   F.stop(t,u) ~ F.stop(T, const Data&)
 */
template<int NDIM, class T = double> class ODErkck;

/************************************************************/
template<int NDIM, class T>
class ODErkck {
public:
  typedef Array1D<T,NDIM> Data;
  struct Soln {
    std::vector<T> t;
    std::vector<Data> u;
    void clear() { t.clear(); u.clear(); }
  };

  /* 
   * return code: 0 success; 1 user stop; 2 minimum step size 
   */
  template<class F>
  int solve(F& _f, T _tmin, T _tmax, const Data& _u0, 
	    T _abstol, T _hmin, T _hmax, Soln& _s) const;
};/* ODErkck */

/************************************************************/
template<int NDIM, class T>
template<class F>
int ODErkck<NDIM,T>::solve(F& _f, T _tmin, T _tmax, const Data& _u0, 
			   T _abstol, T _hmin, T _hmax, Soln& _s) const {
  /*** Define Cash-Karp scheme ***/
  const T c[]   = {0.0, 1.0/5, 3.0/10, 3.0/5, 1.0, 7.0/8};
  const T a0[]  = {0.0};
  const T a1[]  = {1.0/5};
  const T a2[]  = {3.0/40, 9.0/40};
  const T a3[]  = {3.0/10, -9.0/10, 6.0/5};
  const T a4[]  = {-11.0/54, 5.0/2, -70.0/27, 35.0/27};
  const T a5[]  = {1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096};
  const T bO4[] = {2825.0/27648, 0.0, 18575.0/48384, 13525.0/55296, 277.0/14336, 1.0/4};
  const T bO5[] = {37.0/378, 0.0, 250.0/621, 125.0/594, 0.0, 512.0/1771};
  const T* a[6] = {a0,a1,a2,a3,a4,a5};

  /*** Initialize ***/
  T hnew, h = _hmax;
  T tn = _tmin;
  Data un = _u0;
  Data ks[6];
  Data utemp, uO4, uO5;
  T tau; //LTE
  int i, j;

  _s.clear();
  _s.t.push_back(tn);
  _s.u.push_back(un);

  /*** rkck ODE Solver ***/
  while(tn+h < _tmax) {
    //Cash-Karp
    for(i = 0; i < 6; i++) {
      utemp = un;
      for(j = 0; j < i; j++) utemp += a[i][j]*ks[j];
      _f.d(tn+h*c[i], utemp, ks[i]); ks[i] *= h;
    }
    uO4 = un; uO5 = un;
    for(i = 0; i < 6; i++) {
      uO4 += bO4[i]*ks[i];
      uO5 += bO5[i]*ks[i];
    }
    //LTE
    tau = 0.0;
    for(i = 0; i < NDIM; i++) tau = max(tau,abs(uO5[i]-uO4[i]));
    if(tau <= _abstol*h) {//accept the step
      tn += h; un = uO4;
      _s.t.push_back(tn);
      _s.u.push_back(uO5);
      if(_f.stop(tn,un)) return 1;
    }
    //step size control
    hnew = std::pow((0.5*_abstol/tau*std::pow(h,5)),0.25);
    h = min(min(max(hnew,0.1*h),4.0*h),_hmax);
    if(h < _hmin) return 2;
    //adjust for final time
    if(tn+h > _tmax) h = _tmax-tn;
  }//end while(...
  return 0;
}

}/* Util */

#endif
