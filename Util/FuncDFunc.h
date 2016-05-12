#ifndef _UTIL_FUNC_DFUNC_H
#define _UTIL_FUNC_DFUNC_H

#include "util.h"
#include <cstddef>

namespace Util {

template<class T = double> class FuncDFuncBasic;

template<class U, class T = double> class FuncDFunc;

/************************************************************/
template<class T>
class FuncDFuncBasic {
public:
  virtual ~FuncDFuncBasic() {};
  virtual T operator()(T* x) = 0;
  virtual void operator()(T* x, T* g) = 0;
  virtual void operator()(T* x, T* g, T f) = 0;
};

/************************************************************/
template<class U, class T>
class FuncDFunc: public FuncDFuncBasic<T> {
public:
  typedef T (U::*FuncPtr)(T*);
  typedef void (U::*DFuncPtr)(T*,T*);
  typedef void (U::*CallBackPtr)(T*,T*,T);

  FuncDFunc();
  FuncDFunc(U& _u, 
	    int _nDOF,
	    FuncPtr _funcPtr, 
	    DFuncPtr _dfuncPtr, 
	    CallBackPtr _callBackPtr = (CallBackPtr)NULL,
	    std::ostream* _callBackStreamPtr = (std::ostream*)NULL);
  FuncDFunc(const FuncDFunc& fdf);
  FuncDFunc& operator= (const FuncDFunc& fdf);
  T operator()(T* x);
  void operator()(T* x, T* g);
  void operator()(T* x, T* g, T f);
  void init(U& _u, 
	    int _nDOF,
	    FuncPtr _funcPtr, 
	    DFuncPtr _dfuncPtr, 
	    CallBackPtr _callBackPtr = (CallBackPtr)NULL, 
	    std::ostream* _callBackStreamPtr = (std::ostream*)NULL);

protected:
  U* pu;
  int nDOF;
  FuncPtr funcPtr;
  DFuncPtr dfuncPtr;
  CallBackPtr callBackPtr;
  std::ostream* callBackStreamPtr;
  int iter;
};

template<class U, class T>
inline FuncDFunc<U, T>::FuncDFunc() {}

template<class U, class T>
inline FuncDFunc<U, T>::FuncDFunc(U& _u,
				  int _nDOF,
				  FuncPtr _funcPtr, 
				  DFuncPtr _dfuncPtr, 
				  CallBackPtr _callBackPtr, 
				  std::ostream* _callBackStreamPtr) {
  init(_u, _nDOF, _funcPtr, _dfuncPtr, _callBackPtr, _callBackStreamPtr);
}

template<class U, class T>
inline FuncDFunc<U, T>::FuncDFunc(const FuncDFunc& fdf) { *this = fdf; }

template<class U, class T>
FuncDFunc<U, T>& FuncDFunc<U, T>::operator= (const FuncDFunc<U, T>& fdf) {
  if(this != &fdf) {
    pu = fdf.pu;
    nDOF = fdf.nDOF;
    funcPtr = fdf.funcPtr;
    dfuncPtr = fdf.dfuncPtr;
    callBackPtr = fdf.callBackPtr;
    callBackStreamPtr = fdf.callBackStreamPtr;
    iter = fdf.iter;
  }
  return *this;
}

template<class U, class T>
inline T FuncDFunc<U, T>::operator()(T* x) { return (pu->*funcPtr)(x); }

template<class U, class T>
inline void FuncDFunc<U, T>::operator()(T* x, T* g) { (pu->*dfuncPtr)(x,g); }

template<class U, class T>
inline void FuncDFunc<U, T>::operator()(T* x, T* g, T f) {
  if(callBackPtr != (CallBackPtr)NULL) {
    (pu->*callBackPtr)(x,g,f);
  }
  if(callBackStreamPtr != (std::ostream*)NULL) {
    (*callBackStreamPtr) << "<FuncDFunc::callBack> " 
			 << iter++ << " "
			 << f << " " 
			 << norm1(g,nDOF) 
			 << std::endl;
  }
}

template<class U, class T>
void FuncDFunc<U, T>::init(U& _u,
			   int _nDOF,
			   FuncPtr _funcPtr, 
			   DFuncPtr _dfuncPtr, 
			   CallBackPtr _callBackPtr, 
			   std::ostream* _callBackStreamPtr) {
  pu = &_u;
  nDOF = _nDOF;
  funcPtr = _funcPtr;
  dfuncPtr = _dfuncPtr;
  callBackPtr = _callBackPtr;
  callBackStreamPtr = _callBackStreamPtr;
  iter = 0;
} 

}/* Util */

#endif
