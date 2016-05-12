#include "../AMod/Molecule.h"
#include "../Util/RNG.h"
#include <string>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 4) {
    cout << "Usage: testRandRmAtoms <in_mol> <out_mol> <n_rm_atoms>" << endl;
    return 1;
  }  

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  ifstream fin;
  ofstream fout;
  int nRmAtoms;

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

  Util::readWord(argv[3], nRmAtoms);
  if(nRmAtoms < 0 || nRmAtoms > mol.natoms()) {
    cout << "Error: <n_rm_atoms> must be in [0," 
	 << mol.natoms() << "]" << endl;
    return 5;
  }

  /********** REMOVE ATOMS **********/
  Util::RNG::seed(time(NULL));
  int i, j, fixed, atomID;
  for(i = 0; i < nRmAtoms; i++) {
    while(true) {
      atomID = int(Util::RNG::uniform_SFMT()*mol.natoms());
      fixed = 0; for(j = 0; j < 3; j++) fixed += mol.atom(atomID)->fixed[j];
      if(fixed == 0) break;
    }
    mol.deleteAtom(atomID);
  }
  mol.topo().update(AMod::CASE_GENERAL);

  /********** OUTPUT MOL **********/
  mol.io().dumpTxt(fout);

  return 0;
}
