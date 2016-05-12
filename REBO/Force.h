#ifndef _REBO_FORCE_H
#define _REBO_FORCE_H

#include "Energy.h"

namespace REBO {

class Force {
public:
  Force(Energy& energy);
  ~Force();
  void computeBF(const AMod::Arrow& ar, double* f);
  void computeAF(const AMod::Atom& atom, double* f);
  void computeF(AMod::Molecule& mol);

private:
  BondData& bondData;

  Force(const Force& tbforce);
  Force& operator= (const Force& tbforce);
  void compute_f0(const AMod::Arrow& ar0, double* f);
  void compute_f1(const AMod::Arrow& ar0, const AMod::Arrow& ar1, double* f);
  void compute_f2(const AMod::Arrow& ar0, const AMod::Arrow& ar1, const AMod::Arrow& ar2, double* f);
};

} /* REBO */

#endif
