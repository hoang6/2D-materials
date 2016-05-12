#ifndef _MD_LJ_H
#define _MD_LJ_H

#include <cmath>

namespace MD {

class LJ {
public:
  LJ(double _e, double _s);
  virtual ~LJ();
  virtual double energy(double _r) const;
  virtual double denergy(double _r) const;

protected:
  double e, s;

  static double ene(double _e, double _u);
  static double dene(double _e, double _u);
};

class LJCutoff: public LJ {
public:
  LJCutoff(double _e, double _s, double _c);
  virtual double energy(double _r) const ;
  virtual double denergy(double _r) const;

private:
  double c, ic;
  double eneOffset;
  double deneOffset;
};

/************************************************************/
inline LJ::LJ(double _e, double _s): e(_e), s(_s) {}

inline LJ::~LJ() {}

inline double LJ::energy(double _r) const {
  return ene(e,s/_r);
}

inline double LJ::denergy(double _r) const {
  return dene(e,s/_r);
}

inline double LJ::ene(double _e, double _u) {
  return 4.0*_e*(std::pow(_u,12)-std::pow(_u,6)); 
}

inline double LJ::dene(double _e, double _u) {
  return 4.0*_e*(12.0*std::pow(_u,11)-6.0*std::pow(_u,5));
}

/************************************************************/
inline LJCutoff::LJCutoff(double _e, double _s, double _c): LJ(_e,_s),c(_c) {
  ic = 1.0/c;
  eneOffset = LJ::ene(e,ic);
  deneOffset = LJ::dene(e,ic);;
}

inline double LJCutoff::energy(double _r) const {
  double u = s/_r;
  if(u < ic) return 0.0;
  return LJ::ene(e,u)+(ic*ic/u-ic)*deneOffset-eneOffset;
}

inline double LJCutoff::denergy(double _r) const {
  double u = s/_r;
  if(u < ic) return 0.0;
  return -LJ::dene(e,u)*u*u/s+deneOffset*(ic*ic/s);
}

}/* MD */

#endif
