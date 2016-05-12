#include "../AMod/Molecule.h"
#include "../Util/constants.h"
#include <iostream>
#include <cmath>

int main() {
  const double amp = 1.0;
  const double wlen = 8.0;

  AMod::Molecule mol;
  mol.io().sourceTxt(std::cin);
  AMod::MolIO::setFormat(std::cout,12);

  AMod::Molecule molA = mol;  
  for(int i = 0; i < molA.natoms(); i++) {
    molA.atom(i)->x[2] += amp*std::sin(molA.atom(i)->x[0]*2.0*Util::PI/wlen);
  }
  molA.topo().update(AMod::CASE_GENERAL);  
  std::cout << molA;

  AMod::Molecule molB = mol;
  for(int i = 0; i < molB.natoms(); i++) {
    molB.atom(i)->x[2] += amp*std::sin(molB.atom(i)->x[1]*2.0*Util::PI/wlen);
  }
  molB.topo().update(AMod::CASE_GENERAL);  
  std::cout << molB;

  return 0;
}
