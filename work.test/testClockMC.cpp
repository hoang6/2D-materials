#include "../AMod/AMod.h"
#include "../MC/MC.h"
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <vector>

struct Action {
  int atomID;
  double dx[3];
  bool accept;
};

void readxyz(std::istream& _in, AMod::Molecule& _mol);

int main(int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  if(argc != 3) {
    cout << "Usage: testClockMC <file_xyz> <file_MC_action>" << endl;
    return 1;
  }
  
  AMod::Molecule mol;
  std::ifstream mfin(argv[1]);
  if(!mfin) {
    cout << "Error: cannot read file: " << argv[1] << endl;
    return 2;
  }
  readxyz(mfin,mol);
  mfin.close();
  std::ifstream afin(argv[2]);
  if(!afin) {
    cout << "Error: cannot read file: " << argv[2] << endl;
    return 3;
  }
  std::vector<Action> actions;
  while(!afin.eof()) {
    Action a;
    afin >> a.atomID >> a.dx[0] >> a.dx[1] >> a.dx[2] >> a.accept;
    if(!afin.fail()) actions.push_back(a);
  }
  
  MC::MCBasic* pmc = MC::allocMC("LCBOPII");
  pmc->attach(mol);
  pmc->init();
  double new_x[3];
  std::clock_t tic,toc;
  typedef std::vector<Action>::size_type sz_t;
  tic = std::clock();
  for(sz_t kStep = 0; kStep < actions.size(); kStep++) {
    const Action& a = actions[kStep];
    const double* x = mol.atom(a.atomID)->x;
    for(int i = 0; i < 3; i++) 
      new_x[i] = x[i]+a.dx[i];
    pmc->dEnergyMoveAtom(a.atomID,new_x);
    if(!a.accept) pmc->unMoveAtom();
  }
  toc = std::clock();
  pmc->final();
  delete pmc;

  cout.precision(12);
  cout << "Time = " 
       << (double(toc-tic)/CLOCKS_PER_SEC) << " seconds" << endl;
  cout << "Pot  = " << mol.potential() << " eV" << endl;

  /*
  AMod::KernelLCBOPIIN kernel(mol);
  kernel.init();
  kernel.update();
  kernel.energy();
  kernel.final();
  cout << "Ref  = " << mol.potential() << " eV" << endl;
  */

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
