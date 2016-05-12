#ifndef _UTIL_ARRAY1D_H
#define _UTIL_ARRAY1D_H

namespace Util {

template<class T, int N>
class Array1D {
public:
  T data[N];

  Array1D();
  Array1D(const T& _t);
  Array1D(const T _data[N]);
  Array1D(const Array1D& array1D);
  Array1D& operator= (const Array1D& array1D);
  /********** operators **********/
  T& operator[] (int n) { return data[n]; }
  const T& operator[] (int n) const { return data[n]; }
  Array1D<T,N>& operator+= (const Array1D<T,N>& _rhs) {
    for(int i = 0; i < N; i++) data[i] += _rhs.data[i]; return *this; }
  Array1D<T,N>& operator-= (const Array1D<T,N>& _rhs) {
    for(int i = 0; i < N; i++) data[i] -= _rhs.data[i]; return *this; }
  Array1D<T,N>& operator*= (const T& _k) {
    for(int i = 0; i < N; i++) data[i] *= _k; return *this; }
  friend const Array1D<T,N> operator+ (const Array1D<T,N>& _a1, const Array1D<T,N>& _a2) {
    return (Array1D<T,N>(_a1)+=_a2); }
  friend const Array1D<T,N> operator- (const Array1D<T,N>& _a1, const Array1D<T,N>& _a2) {
    return (Array1D<T,N>(_a1)-=_a2); }
  friend const Array1D<T,N> operator* (const T& _k, const Array1D<T,N>& _a) {
    return (Array1D<T,N>(_a)*=_k);
  }
};

/************************************************************/
template<class T, int N>
inline Array1D<T, N>::Array1D() { 
  for(int k = 0; k < N; k++) data[k] = T(0); 
}

template<class T, int N>
inline Array1D<T, N>::Array1D(const T& _t) { 
  for(int k = 0; k < N; k++) data[k] = _t; 
}

template<class T, int N>
inline Array1D<T, N>::Array1D(const T _data[N]) { 
  for(int k = 0; k < N; k++) data[k] = _data[k]; 
}

template<class T, int N>
inline Array1D<T, N>::Array1D(const Array1D<T, N>& array1D) { *this = array1D; }

template<class T, int N>
Array1D<T, N>& Array1D<T, N>::operator= (const Array1D<T, N>& array1D) {
  if(this != &array1D) 
    for(int k = 0; k < N; k++) data[k] = array1D.data[k];
  return *this;
}

} /* Util */

#endif
