#ifndef _UTIL_QUENCH_VERLET_H
#define _UTIL_QUENCH_VERLET_H

#include "Verlet.h"

namespace Util {

template<class F, class T = double> class QuenchVerlet;

template<class F, class T>
bool quenchVerlet(T x[], int n, T tol, int maxIterNum, T dt, 
		  const T mass[], F& funcd, 
		  int& iterNum, int& keyIterNum);

/************************************************************/
template<class F, class T>
class QuenchVerlet: public Verlet<F, T> {
public:
  QuenchVerlet();
  QuenchVerlet(T _tol, int _maxIterNum, T _dt, const T* _mass);
  QuenchVerlet(const QuenchVerlet& quenchVerlet);
  virtual ~QuenchVerlet();
  QuenchVerlet& operator= (const QuenchVerlet& quenchVerlet);

protected:
  virtual void engine();
};

/************************************************************/
template<class F, class T>
inline QuenchVerlet<F, T>::QuenchVerlet() {}

template<class F, class T>
inline QuenchVerlet<F, T>::QuenchVerlet(T _tol, int _maxIterNum, T _dt, const T* _mass): Verlet<F, T>(_tol, _maxIterNum, _dt, _mass) {}

template<class F, class T>
QuenchVerlet<F, T>::QuenchVerlet(const QuenchVerlet<F, T>& quenchVerlet) { *this = quenchVerlet; }

template<class F, class T>
inline QuenchVerlet<F, T>::~QuenchVerlet() {}

template<class F, class T>
QuenchVerlet<F, T>& QuenchVerlet<F, T>::operator= (const QuenchVerlet<F, T>& quenchVerlet) {
  if(this != &quenchVerlet) {
    this->Verlet<F, T>::operator= (quenchVerlet);
  }
  return *this;
}

template<class F, class T>
void QuenchVerlet<F, T>::engine() {
  //quench by projecting v to a
  T anorm2 = vdotv(this->a, this->a, this->n);
  T adotv = vdotv(this->a, this->v, this->n);
  if(adotv > 0)
    vcopy(adotv/anorm2, this->a, this->v, this->n);
  else
    vassign(0.0, this->v, this->n);
}

/************************************************************/
template<class F, class T>
bool quenchVerlet(T x[], int n, T tol, int maxIterNum, T dt, 
		  const T mass[], F& funcd, 
		  int& iterNum, int& keyIterNum) {
  QuenchVerlet<F, T> _verlet(tol, maxIterNum, dt, mass);
  bool flag = _verlet.minimize(x, n, funcd);
  iterNum = _verlet.iterNumber();
  keyIterNum = _verlet.keyIterNumber();
  return flag;
}

}/* Util */

#endif
