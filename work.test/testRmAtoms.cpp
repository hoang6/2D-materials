#include "../AMod/Molecule.h"
#include <string>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc < 4) {
    cout << "Usage: testRmAtoms <in_mol> <out_mol> <rm_atom_IDs>" << endl;
    return 1;
  }  

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  ifstream fin;
  ofstream fout;
  std::vector<int> rmAtomIDs;

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

  for(int i = 3; i < argc; i++) {
    int id;
    Util::readWord(argv[i], id);
    if(id < 0 || id >= mol.natoms()) {
      cout << "Error: atom_ID must be in [0," << mol.natoms() << ")" << endl;
      return 5;
    }
    rmAtomIDs.push_back(id);
  }

  /********** REMOVE ATOMS **********/
  int n = rmAtomIDs.size();
  std::vector<AMod::Atom*> rmAtomPtrs;
  rmAtomPtrs.resize(rmAtomIDs.size());
  for(int i = 0; i < n; i++) rmAtomPtrs[i] = &mol.atom(rmAtomIDs[i]);
  for(int i = 0; i < n; i++)
    mol.deleteAtom(rmAtomPtrs[i]->id());
  mol.topo().update(AMod::CASE_GENERAL);

  /********** OUTPUT MOL **********/
  mol.io().dumpTxt(fout);

  return 0;
}
