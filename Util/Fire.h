#ifndef _UTIL_FIRE_H
#define _UTIL_FIRE_H

#include "Verlet.h"
#include <cmath>

namespace Util {

template<class F, class T = double> class Fire;

template<class F, class T>
bool fire(T x[], int n, T tol, int maxIterNum, T dt, 
	  const T mass[], F& funcd, 
	  int& iterNum, int& keyIterNum);

/************************************************************/
template<class F, class T>
class Fire: public Verlet<F, T> {
public:
  Fire();
  Fire(T _tol, int _maxIterNum, T _dt, const T* _mass);
  Fire(const Fire& fire);
  virtual ~Fire();
  Fire& operator= (const Fire& fire);

protected:
  virtual void initialize(T* _x, int _n, F& _funcd);
  virtual void engine();

  T fInc, fDec, dtMax, aStart, fa, alpha;
  int iterMin, cut;
};

/************************************************************/
template<class F, class T>
inline Fire<F, T>::Fire() {}

template<class F, class T>
inline Fire<F, T>::Fire(T _tol, int _maxIterNum, T _dt, const T* _mass): Verlet<F, T>(_tol, _maxIterNum, _dt, _mass) {}

template<class F, class T>
Fire<F, T>::Fire(const Fire<F, T>& fire) { *this = fire; }

template<class F, class T>
inline Fire<F, T>::~Fire() {}

template<class F, class T>
Fire<F, T>& Fire<F, T>::operator= (const Fire<F, T>& fire) {
  if(this != &fire) {
    this->Verlet<F, T>::operator= (fire);
  }
  return *this;
}

template<class F, class T>
void Fire<F, T>::initialize(T* _x, int _n, F& _funcd) {
  Verlet<F, T>::initialize(_x, _n, _funcd);
  fInc = 1.1;
  fDec = 0.5;
  dtMax = this->deltaTime();
  aStart = 0.1;
  fa = 0.99;
  iterMin = 5;
  alpha = aStart;
  cut = 0;
}

template<class F, class T>
void Fire<F, T>::engine() {
  //FIRE: fast inertial relaxation engine
  T anorm2 = vdotv(this->a, this->a, this->n);
  T vnorm2 = vdotv(this->v, this->v, this->n);
  T adotv = vdotv(this->a, this->v, this->n);
  kdotvadd(1-alpha, this->v, alpha*(std::sqrt(vnorm2)/std::sqrt(anorm2)), this->a, this->v, this->n);
  if(adotv > 0 && this->iter-cut > iterMin) {
    this->deltaTime(min(this->deltaTime()*fInc, dtMax));
    alpha *= fa;
  }
  else if(adotv <= 0) {
    vassign(0.0, this->v, this->n);
    cut = this->iter;
    this->deltaTime(this->deltaTime()*fDec);
    alpha = aStart;
  }
}

/************************************************************/
template<class F, class T>
bool fire(T x[], int n, T tol, int maxIterNum, T dt, 
	  const T mass[], F& funcd, 
	  int& iterNum, int& keyIterNum) {
  Fire<F, T> _verlet(tol, maxIterNum, dt, mass);
  bool flag = _verlet.minimize(x, n, funcd);
  iterNum = _verlet.iterNumber();
  keyIterNum = _verlet.keyIterNumber();
  return flag;
}

} /* Util */

#endif
