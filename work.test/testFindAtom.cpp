#include "../AMod/Molecule.h"
#include <string>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 5) {
    cout << "Usage: testFindAtom <mol_file> <x> <y> <z>" << endl;
    return 1;
  }  

  /********** PARSE ARGUMENT **********/
  AMod::Molecule mol;
  ifstream fin;
  double x[3];

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

  for(int i = 0; i < 3; i++) {
    if(!Util::readWord(argv[i+2], x[i])) {
      cout << "Error: fail to read x y z" << endl;
      return 4;
    }
  }

  /********** Query ATOMS **********/
  int keyid = 0, n = mol.natoms();
  double minDist = 0.0;
  for(int i = 0; i < n; i++) {
    double d = mol.topo().distance(i,x);
    if(i == 0 || d < minDist) {
      keyid = i; minDist = d;
    }
  }
  const double* xa = mol.atom(keyid)->x;
  cout << "atomID = " << keyid << "; x = ("
       << xa[0] << ", " << xa[1] << ", " << xa[2] << "); "
       << "distance = " << minDist << endl;

  return 0;
}
