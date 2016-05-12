#include "../AMod/AMod.h"
#include "../Util/util.h"
#include "../Util/RNG.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 4) {
    cout << "Usage: testOptim <REBO/LCBOPII> <file_mol> <out_file_mol>" << endl;
    return 1;
  }

  /********** Parse Arguments **********/  
  AMod::Molecule mol;
  AMod::MolOptim mopt(mol);
  AMod::PotType potType;
  ifstream fin;
  ofstream fout;

  string potTag = argv[1];
  string fileIn = argv[2];
  string fileOut = argv[3];

  if(potTag == "REBO")
    potType = AMod::REBO;
  else if(potTag == "LCBOPII")
    potType = AMod::LCBOPIIN;
  else {
    cout << "Error: unrecognized potentail: " << potTag << endl;
    return 2;
  }

  fin.open(fileIn.c_str());
  if(!fin.good()) {
    cout << "Error: cannot open file: " << fileIn << endl;
    return 3;
  }
  fin >> mol;
  fin.close();
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 4;
  }

  fout.open(fileOut.c_str());
  if(!fout.good()) {
    cout << "Error: cannot open file: " << fileOut << endl;
    return 5;
  }
  
  /********** Relax **********/
  mopt.init(potType);
  mopt.optimize(AMod::MolOptim::DonlpPara(fileOut,false,0.01));
  mopt.final();

  /********** Output **********/
  AMod::MolO::setFormat(fout);
  fout << mol;  
  fout.close();

  return 0;
}
