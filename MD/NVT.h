#ifndef _MD_NVT_H
#define _MD_NVT_H

#include "NVE.h"

namespace MD {

class NVT: public NVE {//total momentum must = 0! use Data::resetTotalMomentum()
public:
  struct Para;
  NVT();
  virtual ~NVT();
  void reset();
  void init(Data& _data, const Para& _para);
  virtual void step(double _dt);
  int ndof() const { return nf; }
  double eExt() const { return eext; }
  double kTMean() const { return kT; }
  virtual void dumpEnergy(std::ostream& _fout) const;

protected:
  void appNHOper(double _dt);

  double xi, vxi;
  double kT, Q;
  int nf;
  double eext;
};

struct NVT::Para { 
  double kT; double tNH; 
  Para() { reset(); }
  Para(double _kT, double _tNH): kT(_kT), tNH(_tNH) {}
  void reset() { kT = 0.0; tNH = 0.0; }
};

}/* MD */

#endif
