#include "../AMod/Molecule.h"
#include "../AMod/constants.h"
#include "../Util/util.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 7) {
    cout << "Usage: testDispAtom <in_mol> <out_mol> <atom_ID> <disp>" << endl;
    return 1;
  }

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  ifstream fin;
  ofstream fout;
  int atomID;
  double disp[3];

  fin.open(argv[1]);
  if(!fin.good()) {
    cout << "Error: cannot open file: " << argv[1] << endl;
    return 2;
  }
  fin >> mol;
  fin.close();
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 3;
  }
  mol.setStdAxes();

  fout.open(argv[2]);
  if(!fout.good()) {
    cout << "Error: cannot open file: " << argv[2] << endl;
    return 4;
  }
  AMod::MolO::setFormat(fout,12);

  Util::readWord(argv[3], atomID);
  if(atomID < 0 || atomID >= mol.natoms()) {
    cout << "Error: atomID must be in [0," << mol.natoms() << ")" << endl;
    return 5;
  }

  for(int i = 0; i < 3; i++)
    Util::readWord(argv[4+i], disp[i]);

  /********** DISPLACE ATOM **********/
  for(int i = 0; i < 3; i++)
    mol.atom(atomID)->x[i] += disp[i];
  mol.topo().update(AMod::CASE_GENERAL);

  /********** OUTPUT MOL **********/
  mol.io().dumpTxt(fout);

  return 0;
}
