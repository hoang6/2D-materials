#ifndef _UTIL_ARRAY2D_H
#define _UTIL_ARRAY2D_H

namespace Util {

template<class T, int M, int N>
class Array2D {
public:
  T data[M][N];

  Array2D();
  Array2D(const T& _t);
  Array2D(const T& _tDiag, const T& _tOthers);
  Array2D(const T _data[M][N]);
  Array2D(const Array2D& array2D);
  Array2D& operator= (const Array2D& array2D);
  T* operator[] (int m);
  const T* operator[] (int m) const;
  void set_diagonal(const T& _t);
};

/************************************************************/
template<class T, int M, int N>
Array2D<T, M, N>::Array2D() {
  for(int k = 0; k < M; k++)
    for(int j = 0; j < N; j++) data[k][j] = T(0);
}

template<class T, int M, int N>
Array2D<T, M, N>::Array2D(const T& _t) {
  for(int k = 0; k < M; k++)
    for(int j = 0; j < N; j++) data[k][j] = _t;
}

template<class T, int M, int N>
Array2D<T, M, N>::Array2D(const T& _tDiag, const T& _tOthers) {
  for(int k = 0; k < M; k++)
    for(int j = 0; j < N; j++) data[k][j] = (k==j ? _tDiag : _tOthers);
}

template<class T, int M, int N>
Array2D<T, M, N>::Array2D(const T _data[M][N]) {
  for(int k = 0; k < M; k++)
    for(int j = 0; j < N; j++) data[k][j] = _data[k][j];
}

template<class T, int M, int N>
inline Array2D<T, M, N>::Array2D(const Array2D<T, M, N>& array2D) { *this = array2D; }

template<class T, int M, int N>
Array2D<T, M, N>& Array2D<T, M, N>::operator= (const Array2D<T, M, N>& array2D) {
  if(this != &array2D) {
    for(int k = 0; k < M; k++)
      for(int j = 0; j < N; j++) data[k][j] = array2D.data[k][j];
  }
  return *this;
}

template<class T, int M, int N>
inline T* Array2D<T, M, N>::operator[] (int m) { return data[m]; }

template<class T, int M, int N>
inline const T* Array2D<T, M, N>::operator[] (int m) const { return data[m]; }

template<class T, int M, int N>
void Array2D<T, M, N>::set_diagonal(const T& _t) {
  int L = (M<N ? M : N);
  for(int k = 0; k < L; k++) data[k][k] = _t;
}

} /* Util */

#endif
