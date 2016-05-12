#include "NVT.h"
#include "../Util/constants.h"

namespace MD {

NVT::NVT() { reset(); }

NVT::~NVT() {}

void NVT::reset() {
  NVE::reset();
  xi = 0.0;
  vxi = 0.0;
  kT = 0.0;
  Q = 0.0;
  nf = 0;
  eext = 0.0;
}

void NVT::init(Data& _data, const Para& _para) {
  NVE::init(_data);
  xi = 0.0;
  vxi = 0.0;
  kT = _para.kT;
  nf = NVE::ndata()-3;
  Q = nf*kT*(_para.tNH*_para.tNH/(4.0*Util::PI*Util::PI));
}

void NVT::step(double _dt) {
  appNHOper(_dt);
  NVE::step(_dt);
  appNHOper(_dt);
  eext = Q*vxi*vxi/2.0+ndof()*kT*xi;
}

void NVT::appNHOper(double _dt) {
  int i;
  double G, s;

  G = (2.0*NVE::ekin-ndof()*kT)/Q;
  vxi += G*_dt/4.0;
  xi += vxi*_dt/2.0;
  s = std::exp(-vxi*_dt/2.0);
  NVE::ekin *= (s*s);
  for(i = 0; i < NVE::ndata(); i++) vel(i) *= s;
  G = (2.0*NVE::ekin-ndof()*kT)/Q;
  vxi += G*_dt/4.0;
}

void NVT::dumpEnergy(std::ostream& _fout) const {
  _fout << std::showpos << eKin() << "\t" 
	<< ePot() << "\t" << eExt() << std::endl;
}

}/* MD */
