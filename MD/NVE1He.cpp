#include "NVE1He.h"

namespace MD {

NVE1He::NVE1He(): ljc(0.0014,2.74,2.50) { reset(); }

NVE1He::~NVE1He() {}

void NVE1He::reset() {
  NVE::reset();
  topoHeC.reset();
  ehec = 0.0;
  atomHeID = 0;
}

void NVE1He::init(Data& _data, const Para& _para) {
  AMod::Molecule::AData atomHeData(2,_para.posHe);
  atomHeID = _data.mol.newAtom(-1,atomHeData).id();
  for(int i = 0; i < 3; i++) _data.vel.push_back(_para.velHe[i]);
  topoHeC.reset();
  topoHeC.attach(_data.mol);
  topoHeC.sync();
  topoHeC.mode(_data.mol.topo().mode());
  topoHeC.bondRange().insert(2,6,0.00,6.85);
  topoHeC.cellLength() = 6.90;
  topoHeC.update(AMod::CASE_CHANGE_MODE);

  NVE::init(_data);
}

void NVE1He::final() {
  NVE::final();
  topoHeC.reset();
  for(int i = 0; i < 3; i++) (*pvel).pop_back();
  pmol->deleteAtom(atomHeID);
}

void NVE1He::dumpEnergy(std::ostream& _fout) const {
  _fout << std::showpos 
	<< eKin() << "\t" << ePot() << "\t" << eHeC() << std::endl;
}

void NVE1He::computeForce() {
  NVE::computeForce();
  //------------------------------
  int k, kmax = topoHeC.natoms();
  for(k = 0; k < kmax; k++) topoHeC.moveAtom(k);
  AMod::Atom& aHe = topoHeC.atom(atomHeID);
  int i, j;
  double *dx, r, dV, force;
  ehec = 0.0;
  kmax = aHe.arrows.size();
  for(k = 0; k < kmax; k++) {
    dx = aHe.arrows[k].dx;
    r = aHe.arrows[k].bond().r;
    ehec += ljc.energy(r);
    dV = ljc.denergy(r);
    for(j = 0; j < 3; j++) {
      force = -dV/r*(-dx[j]);
      i = 3*aHe.id()+j;
      aHe->fpot[j] += force;
      force = -force;
      i = 3*aHe.arrows[k].moon().id()+j;
      aHe.arrows[k].moon()->fpot[j] += force;
    }
  }
}

}/* MD */
