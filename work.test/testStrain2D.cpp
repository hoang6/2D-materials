#include "../AMod/AMod.h"
#include "../Util/util.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 4) {
    cout << "Usage: testStrain2D <file_mol> <out_file_mol> <Dxx Dyx Dxy Dyy>" << endl;
    return 1;
  }

  /********** Parse Arguments **********/  
  AMod::Molecule mol;
  double Dxx, Dyx, Dxy, Dyy;
  ifstream fin;
  ofstream fout;
  istringstream iss;

  string fileIn = argv[1];
  string fileOut = argv[2];
  string Dtensor = argv[3];

  fin.open(fileIn.c_str());
  if(!fin.good()) {
    cout << "Error: cannot open file: " << fileIn << endl;
    return 2;
  }
  fin >> mol;
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 3;
  }
  fin.close();

  iss.str(Dtensor);
  iss >> Dxx >> Dyx >> Dxy >> Dyy;
  if(iss.fail()) {
    cout << "Error: fail to read D tensor: " << Dtensor << endl;
    return 4;
  }

  fout.open(fileOut.c_str());
  if(!fout.good()) {
    cout << "Error: cannot open file: " << fileOut << endl;
    return 5;
  }

  /********** Relax **********/
  double x, y;
  int j, k, kmax;
  kmax = mol.naxes();
  for(k = 0; k < kmax; k++) {
    x = mol.axis(k).x[0];
    y = mol.axis(k).x[1];
    mol.axis(k).x[0] += x*Dxx+y*Dyx;
    mol.axis(k).x[1] += x*Dxy+y*Dyy;
    for(j = 0; j < 3; j++)
      mol.axis(k).fixed[j] = 1;
  }
  kmax = mol.natoms();
  for(k = 0; k < kmax; k++) {
    x = mol.atom(k)->x[0];
    y = mol.atom(k)->x[1];
    mol.atom(k)->x[0] += x*Dxx+y*Dyx;
    mol.atom(k)->x[1] += x*Dxy+y*Dyy;
  }
  mol.topo().update(AMod::CASE_GENERAL);

  /********** Output **********/
  AMod::MolO::setFormat(fout);
  fout << mol;  
  fout.close();

  return 0;
}
