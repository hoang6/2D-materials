#ifndef _UTIL_VERLET_H
#define _UTIL_VERLET_H

#include "MinimizerBasic.h"
#include "Accessors.h"
#include "util.h"

namespace Util {

template<class F, class T = double> class Verlet;

template<class F, class T>
bool verlet(T x[], int n, T tol, int maxIterNum, T dt, 
	    const T mass[], F& funcd, 
	    int& iterNum, int& keyIterNum);

/************************************************************/
template<class F, class T>
class Verlet: public MinimizerBasic<F, T> {
public:
  Verlet();
  Verlet(T _tol, int _maxIterNum, T _dt, const T* _mass);
  Verlet(const Verlet& verlet);
  virtual ~Verlet();
  Verlet& operator= (const Verlet& verlet);
  virtual void reset();
  virtual bool minimize(T* _x, int _n, F& _funcd);
  virtual typename MinimizerBasic<F, T>::Style style() const;
  int keyIterNumber() const;
  T keyNormGrad() const;

  Accessors<const T*> mass;

protected:
  virtual void initialize(T* _x, int _n, F& _funcd);
  virtual void finalize(int _iter);
  virtual void engine();
  virtual void integral();
  void step();
  void comp_a(T* _x, T* _a);

  T* x; //x(t)
  T* v; //dx/dt
  T* g;
  T* a; //d2x/dt2
  T* massOne;
  int n;
  T t;
  int iter;
  int keyIterNum;
  T keyNormG;
  F* pfuncd;

private:
  void _reset();
};
/************************************************************/
template<class F, class T>
inline Verlet<F, T>::Verlet() { _reset(); }

template<class F, class T>
inline Verlet<F, T>::Verlet(T _tol, int _maxIterNum, T _dt, const T* _mass): MinimizerBasic<F, T>(_tol, _maxIterNum, _dt) { mass(_mass); }

template<class F, class T>
inline Verlet<F, T>::Verlet(const Verlet<F, T>& verlet) { *this = verlet; }

template<class F, class T>
inline Verlet<F, T>::~Verlet() {}

template<class F, class T>
Verlet<F, T>& Verlet<F, T>::operator= (const Verlet<F, T>& verlet) {
  if(this != &verlet) {
    this->MinimizerBasic<F, T>::operator= (verlet);
    mass = verlet.mass;
  }
  return *this;
}

template<class F, class T>
inline void Verlet<F, T>::reset() { _reset(); }

template<class F, class T>
bool Verlet<F, T>::minimize(T* _x, int _n, F& _funcd) {
  T normg(0);
  initialize(_x, _n, _funcd);
  (*pfuncd)(x, g, T(0));
  for(iter = 1; iter <= this->maxIterNumber(); iter++) {
    //adaptive run
    step();    
    (*pfuncd)(x, g, T(0));
    //stop?
    normg = norm1(g, n);
    if(iter == 1 || keyNormG > normg) {
      keyIterNum = iter;
      keyNormG = normg;
    }
    if(this->isbad(normg)) {
      finalize(iter);
      return false;
    }
    if(normg < this->tolerance()) { 
      finalize(iter);
      return true;
    }
  } 

#ifndef NR_SILENT
  std::cerr << "Too many iterations in VV::minimize, "
	    << "|gradient| = " << normg << std::endl;
#endif

  finalize(iter-1);
  return false;
}

template<class F, class T>
inline typename MinimizerBasic<F, T>::Style Verlet<F, T>::style() const { return MinimizerBasic<F, T>::DFUNC_ONLY; }

template<class F, class T>
inline int Verlet<F, T>::keyIterNumber() const { return keyIterNum; }

template<class F, class T>
inline T Verlet<F, T>::keyNormGrad() const { return keyNormG; }

template<class F, class T>
void Verlet<F, T>::initialize(T* _x, int _n, F& _funcd) {
  x = _x;
  n = _n;
  pfuncd = &_funcd;
  t = 0.0;
  v = new T [n];
  g = new T [n];
  a = new T [n];
  massOne = (T*)NULL;
  if(mass() == (T*)NULL) {
    massOne = new T[n];
    for(int k = 0; k < n; k++) massOne[k] = 1.0;
    mass(massOne);
  }
  vassign(0.0, v, n);
  comp_a(x, a);
}

template<class F, class T>
void Verlet<F, T>::finalize(int _iter) {
  if(massOne != (T*)NULL) {
    delete [] massOne;
    mass((T*)NULL);
  }
  delete [] a;
  delete [] g;
  delete [] v;

  this->iterNum = _iter;
}

template<class F, class T>
void Verlet<F, T>::engine() {}

template<class F, class T>
void Verlet<F, T>::integral() {
  for(int k = 0; k < n; k++) {
    v[k] = v[k]+a[k]*this->deltaTime()/2;
    x[k] = x[k]+v[k]*this->deltaTime();
  }
  comp_a(x, a);
  for(int k = 0; k < n; k++)
    v[k] += a[k]*this->deltaTime()/2;
  t += this->deltaTime();
}

template<class F, class T>
void Verlet<F, T>::step() {
  integral();
  engine();
}

template<class F, class T>
void Verlet<F, T>::comp_a(T* _x, T* _a) {
  (*pfuncd)(_x, g);
  for(int k = 0; k < n; k++) _a[k] = -g[k]/mass()[k];
}

template<class F, class T>
void Verlet<F, T>::_reset() {
  this->tolerance(0.05);
  this->maxIterNumber(10000);
  this->deltaTime(0.10);
  this->mass(NULL);
}

/************************************************************/
template<class F, class T>
bool verlet(T x[], int n, T tol, int maxIterNum, T dt, 
	    const T mass[], F& funcd, 
	    int& iterNum, int& keyIterNum) {
  Verlet<F, T> _verlet(tol, maxIterNum, dt, mass);
  bool flag = _verlet.minimize(x, n, funcd);
  iterNum = _verlet.iterNumber();
  keyIterNum = _verlet.keyIterNumber();
  return flag;
}

}/* Util */

#endif
