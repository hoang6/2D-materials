#include "../AMod/Molecule.h"
#include <iostream>

int main() {
  AMod::Molecule mol;
  std::cin >> mol;
  if(std::cin.fail() || mol.natoms() == 0) {
    std::cout << "Error: no molecule is read" << std::endl;
    return 1;
  }
  
  double cx[3];
  const double* lA = mol.axis(0).x;
  const double* lB = mol.axis(1).x;
  cx[0] = (lA[0]+lB[0])/2.0;
  cx[1] = (lA[1]+lB[1])/2.0;
  for(int i = 0; i < mol.natoms(); i++) {
    cx[2] = mol.atom(i)->x[2];
    mol.topo().around(cx,i);
  }

  AMod::MolO::setFormat(std::cout,12);
  std::cout << mol;

  return 0;
}
