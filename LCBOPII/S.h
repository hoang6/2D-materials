#ifndef _LCBOPII_S_H
#define _LCBOPII_S_H

namespace LCBOPII {

double fSd(double x, double p);

double fSu(double x, double p);

double fSdsr(double r);

double fSdlr(double r);

double fSddb(double xdb);

double fSumr(double r);

double fSuM(double N);

double fSuel(double N);

double fSugamma0(double gamma);

double fSdN(double r);

double fSdsat(double N);

double fSugamma2(double gamma);

/************************************************************/
inline double fSd(double x, double p) {
  if(x <= 0.0) return 1.0;
  if(x >= 1.0) return 0.0;
  return (1.0+2.0*x+p*x*x)*(1.0-x)*(1.0-x);
}

inline double fSu(double x, double p) { 
  return 1.0-fSd(x,p);
}

inline double fSdsr(double r) {
  return fSd((r-1.7)/0.5,3.0);
}

inline double fSdlr(double r) {
  return fSd((r-5.5)/0.5,0.0);
}

inline double fSddb(double xdb) {
  return fSd(xdb,0.0);
}

inline double fSumr(double r) {
  return fSu((r-1.7)/0.5,+2.0);
}

inline double fSuM(double N) {
  return fSu((N-2.0),0.0);
}

inline double fSuel(double N) {
  return fSu((N-3.0)/0.2,0.0);
}

inline double fSugamma0(double gamma) {
  return fSu((gamma-0.34)/0.59,0.0);
}

inline double fSdN(double r) {
  return fSd((r-1.7)/0.5,-3.0);
}

inline double fSdsat(double N) {
  return fSd((N-3.0),0.0);
}

inline double fSugamma2(double gamma) {
  return fSu((gamma-0.30)/0.63,0.0);
}

}/* LCBOPII */

#endif
