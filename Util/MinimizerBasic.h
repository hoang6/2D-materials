#ifndef _UTIL_MINIMIZER_BASIC_H
#define _UTIL_MINIMIZER_BASIC_H

#include "Accessors.h"
#include <string>
#include <iostream>

namespace Util {

/*
  class F {
    T operator() (T* x);
    void operator() (T* x, T* grad);
    void operator() (T* x, T* grad, T f);
  };
*/
template<class F, class T = double> class MinimizerBasic;

/************************************************************/
template<class F, class T>
class MinimizerBasic {
public:
  enum Style { DFUNC_ONLY, FUNC_DFUNC };

  MinimizerBasic();
  MinimizerBasic(T _tol, int _maxIterNum, T _dt);
  MinimizerBasic(const MinimizerBasic& minimizerBasic);
  virtual ~MinimizerBasic();
  MinimizerBasic& operator= (const MinimizerBasic& minimizerBasic);
  virtual void reset() = 0;
  virtual bool minimize(T* _x, int _n, F& _funcd) = 0;
  virtual typename MinimizerBasic<F, T>::Style style() const = 0;
  bool isbad(const T& _x) const;
  int iterNumber() const;
  virtual int sourceTxt(std::istream& fin);
  virtual void dumpTxt(std::ostream& fout) const;

  Accessors<T> tolerance;
  Accessors<int> maxIterNumber;
  Accessors<T> deltaTime;

protected:
  int iterNum;
};

/************************************************************/
template<class F, class T>
inline MinimizerBasic<F, T>::MinimizerBasic() {}

template<class F, class T>
MinimizerBasic<F, T>::MinimizerBasic(T _tol, int _maxIterNum, T _dt) {
  tolerance(_tol);
  maxIterNumber(_maxIterNum);
  deltaTime(_dt);
}

template<class F, class T>
inline MinimizerBasic<F, T>::MinimizerBasic(const MinimizerBasic<F, T>& minimizerBasic) { *this = minimizerBasic; }

template<class F, class T>
inline MinimizerBasic<F, T>::~MinimizerBasic() {}

template<class F, class T>
MinimizerBasic<F, T>& MinimizerBasic<F, T>::operator= (const MinimizerBasic<F, T>& minimizerBasic) {
  if(this != &minimizerBasic) {
    tolerance = minimizerBasic.tolerance;
    maxIterNumber = minimizerBasic.maxIterNumber;
    deltaTime = minimizerBasic.deltaTime;
  }
  return *this;
}

template<class F, class T>
inline bool MinimizerBasic<F, T>::isbad(const T& _x) const {
  return (isnan(_x) || isinf(_x));
}

template<class F, class T>
inline int MinimizerBasic<F, T>::iterNumber() const { return iterNum; }

template<class F, class T>
int MinimizerBasic<F, T>::sourceTxt(std::istream& fin) {
  std::string tword;
  if(!fin.good()) return 1;
  fin >> tword >> tolerance();
  fin >> tword >> maxIterNumber();
  fin >> tword >> deltaTime();
  return (fin.good() ? 0 : 1);
}

template<class F, class T>
void MinimizerBasic<F, T>::dumpTxt(std::ostream& fout) const {
  using namespace std;
  fout << "tolerance " << tolerance() << endl;
  fout << "maxIterNumber " << maxIterNumber() << endl;
  fout << "deltaTime " << deltaTime() << endl;
}

}/* Util */

#endif
