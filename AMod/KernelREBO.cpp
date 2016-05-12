#include "KernelREBO.h"
#include "Molecule.h"

namespace AMod {

KernelREBO::KernelREBO() {}

KernelREBO::KernelREBO(Molecule& _mol): KernelBasic(_mol) {}

KernelREBO::~KernelREBO() {}

void KernelREBO::init() {
  pEnergyCalc = new REBO::Energy;
  pForceCalc = new REBO::Force(*pEnergyCalc);
}

void KernelREBO::final() {
  delete pForceCalc;
  delete pEnergyCalc;
}

void KernelREBO::update() {
  Molecule& mol = molecule();
  mol.topo().update(AMod::CASE_AXES);
  mol.topo().update(AMod::CASE_GENERAL);
}

void KernelREBO::update(const double* _x) {
  AMod::Molecule& mol = molecule();
  int k, kmax = mol.natoms();
  for(k = 0; k < kmax; k++) mol.moveAtom(k,_x+3*k);
}

double KernelREBO::energy() {
  Molecule& mol = molecule();
  mol.potential(pEnergyCalc->computeE(mol.topo()));
  mol.work(0.0);
  return mol.totalEnergy();
}

void KernelREBO::denergy() {
  pForceCalc->computeF(molecule());
}

}/* AMod */
