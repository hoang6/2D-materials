#ifndef _LCBOPII_HANDY_H
#define _LCBOPII_HANDY_H

#include "globalConstants.h"
#include <cmath>

namespace LCBOPII {

double square(double x);

double costh(const double a[3], const double b[3]);

double norm(const double a[3]);

void normalize(const double a[3], double b[3]);

void cross(const double a[3], const double b[3], double c[3]);

template<typename T>
T min(const T& a, const T& b);

template<typename T>
T max(const T& a, const T& b);

/************************************************************/
inline double square(double x) { return x*x; }

inline double costh(const double a[3], const double b[3]) {
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])/(norm(a)*norm(b));
}

inline double norm(const double a[3]) {
  double n = a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
  return std::sqrt(n);
}

inline void normalize(const double a[3], double b[3]) {
  double n = norm(a);
  if(n < Const::EPS_DOUBLE) return;
  for(int k = 0; k < 3; k++) b[k] = a[k]/n;
}

inline void cross(const double a[3], const double b[3], double c[3]) {
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = b[0]*a[2]-b[2]*a[0];
  c[2] = a[0]*b[1]-b[0]*a[1];
}

template<typename T>
inline T min(const T& a, const T& b) {
  return (a<b?a:b);
}

template<typename T>
inline T max(const T& a, const T& b) {
  return (a<b?b:a);
}

}/* LCBOPII */

#endif
