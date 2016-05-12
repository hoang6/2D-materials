#include "../AMod/AMod.h"
#include "../SW/SW.h"
#include "../Util/util.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 4) {
    cout << "Usage: testSW <file_mol> <rot_bond> <out_file_mol>" << endl;
    return 1;
  }
  
  /********** Parse Arguments **********/
  AMod::Molecule mol;
  ifstream fin;
  int bondID;
  SW::SWStimulus swStim;
  ofstream fout;

  fin.open(argv[1]);
  if(!fin.good()) {
    cout << "Error: cannot open file " << argv[1] << endl;
    return 2;
  }
  fin >> mol;
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 3;
  }
  fin.close();
  mol.topo().mode(AMod::MODE_FREEZE);
  mol.topo().update(AMod::CASE_CHANGE_MODE);
  if(!Util::readWord(argv[2], bondID)) {
    cout << "Error: fail to read bondID" << endl;
    return 4;
  }
  swStim.set(bondID);
  if(!swStim.commit(mol)) {
    cout << "Error: illegal bondID: " << bondID << endl;
    return 5;
  }
  fout.open(argv[3]);
  if(!fout.good()) {
    cout << "Error: cannot open file: " << argv[3] << endl;
    return 6;
  }
  
  /********** Rotate Bond **********/
  SW::SWCreator sw(mol);
  sw.act(swStim);
  
  /********** Output **********/
  AMod::MolO::setFormat(fout);
  fout << mol;
  fout.close();

  return 0;
}
