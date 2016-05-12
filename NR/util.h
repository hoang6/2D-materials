#ifndef _NR_UTIL_H
#define _NR_UTIL_H

#include <cmath>
#include <iostream>

namespace NR {

/************************************************************/
template<class T> T MIN(const T& A, const T& B);

template<class T> T MAX(const T& A, const T& B);

template <class T> T ABS(const T& a);

template<class T> T SIGN(const T& a, const T& b);

template<class T> T SQR(const T& a);

template<class T> T* vector(long nl, long nh);

template<class T> T** matrix(long nrl, long nrh, long ncl, long nch);

template<class T> void free_vector(T* v, long nl, long nh);

template<class T> void free_matrix(T** m, long nrl, long nrh, long ncl, long nch);

template<class T>
class Vector {
public:
  Vector(long _nl, long _nh);
  ~Vector();
  T* pointer();
private:
  Vector(const Vector& vector);
  Vector& operator= (const Vector& vector);
  long nl;
  long nh;
  T *v;
};

template<class T>
class Matrix {
public:
  Matrix(long _nrl, long _nrh, long _ncl, long _nch);
  ~Matrix();
  T** pointer();
private:
  Matrix(const Matrix& matrix);
  Matrix& operator= (const Matrix& matrix);
  long nrl;
  long nrh;
  long ncl;
  long nch;
  T **m;
};

/************************************************************/
template<class T>
T MIN(const T& A, const T& B) { return ((A) < (B)) ? (A) : (B); }

template<class T>
T MAX(const T& A, const T& B) { return ((A) > (B)) ? (A) : (B); }

template <class T>
T ABS(const T& a) { return ((a) < T(0) ? (-(a)) : (a)); }

template<class T>
T SIGN(const T& a, const T& b) { return ((b) >= T(0) ? ABS(a) : (-ABS(a))); }

template<class T>
T SQR(const T& a) { return a*a; }

template<class T>
T* vector(long nl, long nh) {
  T* v = new T [nh-nl+1+1];
  return v-nl+1;
}

template<class T>
T** matrix(long nrl, long nrh, long ncl, long nch) {
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  T** m;
  
  /* allocate pointers to rows */
  m = new T* [nrow+1];
  m += 1;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl] = new T [nrow*ncol+1];
  m[nrl] += 1;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	
  /* return pointer to array of pointers to rows */
  return m;
}

template<class T>
void free_vector(T* v, long nl, long nh) {
  delete [] (v+nl-1);
}

template<class T>
void free_matrix(T** m, long nrl, long nrh, long ncl, long nch) {
  delete [] (m[nrl]+ncl-1);
  delete [] (m+nrl-1);
}

template<class T>
Vector<T>::Vector(long _nl, long _nh) { 
  nl = _nl; nh = _nh;
  v = vector<T>(nl,nh);
}

template<class T>
Vector<T>::~Vector() { free_vector<T>(v,nl,nh); }

template<class T>
T* Vector<T>::pointer() { return v; }

template<class T>
Matrix<T>::Matrix(long _nrl, long _nrh, long _ncl, long _nch) {
  nrl = _nrl; nrh = _nrh; ncl = _ncl; nch = _nch;
  m = matrix<T>(nrl,nrh,ncl,nch);
}

template<class T>
Matrix<T>::~Matrix() { free_matrix<T>(m,nrl,nrh,ncl,nch); }

template<class T>
T** Matrix<T>::pointer() { return m; }

}/* NR */

#endif
