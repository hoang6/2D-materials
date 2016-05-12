#include "Data.h"

namespace MD {

void Data::resetTotalMomentum() {
  forceConsistency();
  int i, j, natoms = mol.natoms();
  double m, totMass, mv[3];
  totMass = 0.0;
  for(j = 0; j < 3; j++) mv[j] = 0.0;
  for(i = 0; i < natoms; i++) { 
    m = mass(i);
    totMass += m;
    for(j = 0; j < 3; j++) mv[j] += m*vel[3*i+j];
  }
  for(j = 0; j < 3; j++) mv[j] /= totMass;
  for(i = 0; i < natoms; i++)
    for(j = 0; j < 3; j++) vel[3*i+j] -= mv[j];
}

void Data::rescaleVel(double _kT) {
  forceConsistency();
  int i, imax = vel.size();
  double ekin = 0.0; //Ek_tot = 3/2*Na*kT
  for(i = 0; i < imax; i++) ekin += 0.5*mass(i/3)*vel[i]*vel[i];
  double scale = std::sqrt(3.0/2.0*mol.natoms()*_kT/ekin);
  for(i = 0; i < imax; i++) vel[i] *= scale;
}

void Data::rmFreeC(const std::vector<int>& _freeAtomIDs) {
  forceConsistency();
  int i, j, m, n;

  Util::PtrCont<Point> xvel;
  n = mol.natoms();
  xvel.resize(n);
  for(i = 0; i < n; i++)
    for(j = 0; j < 3; j++) xvel[i].x[j] = vel[3*i+j];

  std::vector<AMod::Atom*> atomPtrs;
  n = _freeAtomIDs.size(); atomPtrs.resize(n);
  for(i = 0; i < n; i++) 
    atomPtrs[i] = &mol.atom(i);
  for(i = 0; i < n; i++) {
    m = atomPtrs[i]->id();
    mol.deleteAtom(m);
    xvel.erase(m);
  }

  n = mol.natoms();
  vel.resize(3*n);
  for(i = 0; i < n; i++)
    for(j = 0; j < 3; j++) vel[3*i+j] = xvel[i].x[j];
}

}/* MD */
