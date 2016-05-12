#include "../AMod/AMod.h"
#include <iostream>
#include <string>
#include <ctime>

#define N_REPEAT 200

void readxyz(std::istream& _in, AMod::Molecule& _mol);

int main() {
  AMod::Molecule mol;
  readxyz(std::cin,mol);

  std::clock_t tic,toc;
  AMod::KernelLCBOPIIN *pk = new AMod::KernelLCBOPIIN(mol);
  pk->init();
  tic = std::clock();
  for(int k = 0; k < N_REPEAT; k++) pk->energy();
  toc = std::clock(); 
  int natoms = mol.natoms();
  double esr = pk->kernelF().kernelE().getESR();
  double emr = pk->kernelF().kernelE().getEMR();
  double elr = pk->kernelF().kernelE().getELR();
  double etot = pk->kernelF().kernelE().getEtot();
  pk->final();
  delete pk;

  using std::cout;
  using std::endl;
  cout.precision(16);
  cout << "np = " << natoms << endl;
  cout << "etot, etotpp = " << etot << " " << etot/natoms << endl;
  cout << "esr, esrpp   = " << esr << " " << esr/natoms << endl;
  cout << "emr, emrpp   = " << emr << " " << emr/natoms << endl;
  cout << "elr, elrpp   = " << elr << " " << elr/natoms << endl;
  cout << "time         = " 
       << (double(toc-tic)/CLOCKS_PER_SEC)/N_REPEAT << " seconds" << endl;

  return 0;
}

void readxyz(std::istream& _in, AMod::Molecule& _mol) {
  int i, k, natoms;
  double x[3];
  std::string word;
  
  _mol.reset();
  _mol.topo().mode(AMod::MODE_CELL);
  _mol.topo().bondRange() = LCBOPII::Energy::bondRangeSR();
  _mol.topo().cellLength() = LCBOPII::Energy::cellLengthSR();
  
  _in >> natoms;
  _mol.axes().resize(3);
  _mol.atomsData().resize(natoms);

  for(i = 0; i < 3; i++) _in >> x[i];
#ifdef DOUBLE_EDGE
  for(i = 0; i < 3; i++) x[i] *= 2.0;
#endif
  _mol.axis(0) = AMod::Atom::Data(0,x[0],0.0,0.0);
  _mol.axis(1) = AMod::Atom::Data(0,0.0,x[1],0.0);
  _mol.axis(2) = AMod::Atom::Data(0,0.0,0.0,x[2]);
  for(k = 0; k < natoms; k++) {
    _in >> word;
    for(i = 0; i < 3; i++) _in >> x[i];
    _mol.atomsData()[k] = AMod::Atom::Data(6,x[0],x[1],x[2]);
  }
  _mol.topo().sync();
  _mol.topo().update(AMod::CASE_AXES);
  _mol.topo().update(AMod::CASE_BRUTAL);
}
