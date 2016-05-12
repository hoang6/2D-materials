#include "Energy.h"
#include "S.h"
#include "globalConstants.h"
#include <cmath>

namespace LCBOPII {
  
Energy::Energy() {}

Energy::Energy(const Energy& energy) { *this = energy; }
 
Energy::~Energy() {}

Energy& Energy::operator= (const Energy& energy) {
  if(this != &energy) {
    bondDataSR = energy.bondDataSR;
    bondDataLR = energy.bondDataLR;
    atomDataMR = energy.atomDataMR;
  }
  return *this;
}

void Energy::reset() {
  bondDataSR.clear();
  bondDataLR.clear();
  atomDataMR.clear();
  potSR = 0.0;
  potMR = 0.0;
  potLR = 0.0;
}

AMod::BondRange Energy::bondRangeSR() {
  AMod::BondRange br;
  br.insert(6, 6, 0.0, 2.2);
  return br;
}

AMod::BondRange Energy::bondRangeMR() {
  AMod::BondRange br;
  br.insert(6, 6, 1.7, 4.0);
  return br;
}

AMod::BondRange Energy::bondRangeLR() {
  AMod::BondRange br;
  br.insert(6, 6, 1.7, 6.0);
  return br;
}

double Energy::cellLengthSR() { return 2.25; }

double Energy::cellLengthMR() { return 4.05; }

double Energy::cellLengthLR() { return 6.05; }

void Energy::compute_r0(const AMod::Bond& bondSR) {//SdN
  BondDatumSR& bdSR = bondDataSR[bondSR.id()];
  bdSR.SdN = fSdN(bondSR.r);
  bdSR.fullNeighb = (bdSR.SdN == 1.0 ? 1 : 0);
}

void Energy::compute_r1(const AMod::Atom& atomSR) {//N_
  double N = 0.0;
  int k, kmax = atomSR.arrows.size();
  for(k = 0; k < kmax; k++) 
    N += bondDataSR[atomSR.arrows[k].bond().id()].SdN;
  for(k = 0; k < kmax; k++) {
    BondDatumSR& bdSR = bondDataSR[atomSR.arrows[k].bond().id()];
    int pos = atomSR.arrows[k].position();
    bdSR.N_[pos] = N-bdSR.SdN;
    bdSR.M_[pos] = fSuM(bdSR.N_[pos]);
    bdSR.Suel_[pos] = fSuel(bdSR.N_[pos]);
  }
}

void Energy::compute_r2(const AMod::Bond& bondSR) {//SVsr
  bondDataSR[bondSR.id()].E = fSVsr(bondSR);
}

void Energy::compute_r3(const AMod::Bond& bondLR) {//SVlr
  bondDataLR[bondLR.id()].E = fSVlr(bondLR);
}

void Energy::compute_r4(const AMod::Atom& atomMR,
			const AMod::MolTopo& topoSR) {//Emr
  double sum, v;
  int k, kmax = atomMR.arrows.size();
  sum = 0.0;
  for(k = 0; k < kmax; k++) {
    const AMod::Arrow& ar = atomMR.arrows[k];
    v = fSVmr(ar,topoSR);
    sum += v*v;
  }
  atomDataMR[atomMR.id()].E = -0.5*std::sqrt(sum);
}

double Energy::computeE(const AMod::MolTopo& topoSR,
			const AMod::MolTopo& topoMR,
			const AMod::MolTopo& topoLR) {
  int k;
  int natomsSR = topoSR.natoms();
  int natomsMR = topoMR.natoms();
  int nbondsSR = topoSR.nbonds();
  int nbondsLR = topoLR.nbonds();

  bondDataSR.resize(nbondsSR);
  bondDataLR.resize(nbondsLR);
  atomDataMR.resize(natomsMR);
  
  for(k = 0; k < nbondsSR; k++) compute_r0(topoSR.bond(k));
  for(k = 0; k < natomsSR; k++) compute_r1(topoSR.atom(k));
  for(k = 0; k < nbondsSR; k++) compute_r2(topoSR.bond(k));
  for(k = 0; k < nbondsLR; k++) compute_r3(topoLR.bond(k));
  for(k = 0; k < natomsMR; k++) compute_r4(topoMR.atom(k),topoSR);
  for(k = 0, potSR = 0.0; k < nbondsSR; k++) potSR += bondDataSR[k].E;
  for(k = 0, potLR = 0.0; k < nbondsLR; k++) potLR += bondDataLR[k].E;
  for(k = 0, potMR = 0.0; k < natomsMR; k++) potMR += atomDataMR[k].E;

  return potSR+potLR+potMR;
}

}/* LCBOPII */
