#include "../AMod/Molecule.h"
#include "../AMod/constants.h"
#include "../Util/util.h"
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace std;

void makeGPHnLa(AMod::Molecule& _mol, int _nlayers, double _dist);

int main(int argc, char* argv[]) {
  if(argc != 4) {
    cout << "Usage: testGPHnLa <in_mol> <out_mol> <n_layers>" << endl;
    cout << "       input GPH should be 'armchair'-type" << endl;
    return 1;
  }

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  ifstream fin;
  ofstream fout;
  int nlayers;

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

  Util::readWord(argv[3], nlayers);
  if(nlayers < 1) {
    cout << "Error: n_layers = " << nlayers << " < 1" << endl;
    return 5;
  }

  /********** MAKE DOUBLE GPH **********/
  makeGPHnLa(mol,nlayers,3.35);
  mol.io().dumpTxt(fout);

  return 0;
}

void makeGPHnLa(AMod::Molecule& _mol, int _nlayers, double _dist) {
  AMod::Atom::Data tdata;
  int c, k, j, n = _mol.natoms();
  for(c = 1; c < _nlayers; c++) {
    for(k = 0; k < n; k++) {
      tdata = _mol.atom(k).data();
      tdata.x[0] -= (c%2)*AMod::CC_BOND_LENGTH;
      tdata.x[2] -= c*_dist;
      if(tdata.x[0] < 0) {
	for(j = 0; j < 3; j++)
	  tdata.x[j] += _mol.axis(0).x[j];
      }
      _mol.newAtom(-1, tdata);
    }
  }
}
