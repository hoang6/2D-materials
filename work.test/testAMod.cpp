#include "../AMod/Molecule.h"
#include "../AMod/MolChecker.h"
#include "../AMod/MolAdjuster.h"
#include "../Util/util.h"
#include "../Util/RNG.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 2 && argc != 3) {
    cout << "Usage: testAMod <file_mol> <seed>" << endl;
    return 1;
  }

  AMod::Molecule mol;
  AMod::MolChecker mchecker(mol);
  AMod::MolAdjuster madjuster(mol);
  AMod::Atom::Data data;
  bool flag;
  int i, j, n, m, seed;
  ifstream fin;

  fin.open(argv[1]);
  if(!fin.good()) {
    cout << "Error: cannot open file: " << argv[1] << endl;
    return 2;
  }
  fin >> mol;
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 3;
  }
  fin.close();
  if(argc == 3) {
    if(!Util::readWord(argv[2], seed))
      cout << "Error: seed error" << endl;
    return 4;
  }
  else
    seed = time(NULL);

  Util::RNG::seed(seed);
  cout << "seed = " << seed << endl;
  cout << endl;

  cout << "*********************************************" << endl;
  cout << "Check 1: free" << endl;
  cout << "*********************************************" << endl;
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 4;

  cout << "*********************************************" << endl;
  cout << "Check 2: change topo mode (null)" << endl;
  cout << "*********************************************" << endl;
  mol.topo().mode(AMod::MODE_NULL);
  mol.topo().update(AMod::CASE_CHANGE_MODE);
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 5;

  cout << "*********************************************" << endl;
  cout << "Check 3: del # atoms" << endl;
  cout << "*********************************************" << endl;
  m = mol.natoms();
  n = m*(1.0-Util::RNG::uniform_std());
  cout << "# = " << n << endl;
  for(i = 0; i < n; i++) {
    mol.deleteAtom(mol.natoms()*Util::RNG::uniform_std());
  }
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 6;

  cout << "*********************************************" << endl;
  cout << "Check 4: new # atoms" << endl;
  cout << "*********************************************" << endl;
  n = m*(1.0-Util::RNG::uniform_std());
  cout << "# = " << n << endl;
  for(i = 0; i < n; i++) {
    data.type = AMod::PTE::Carbon;
    for(j = 0; j < 3; j++) data.x[j] = Util::RNG::uniform_std()*20.0;
    mol.newAtom(-1,data);
  }
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 7;

  cout << "*********************************************" << endl;
  cout << "Check 5: move # atoms" << endl;
  cout << "*********************************************" << endl;
  n = mol.natoms()*(1.0-Util::RNG::uniform_std());
  cout << "# = " << n << endl;
  for(i = 0; i < n; i++) {
    for(j = 0; j < 3; j++) data.x[j] = Util::RNG::uniform_std()*20.0;
    mol.moveAtom(mol.natoms()*Util::RNG::uniform_std(),data.x);
  }
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 8;

  cout << "*********************************************" << endl;
  cout << "Check 6: change topo mode (cell)" << endl;
  cout << "*********************************************" << endl;
  mol.topo().mode(AMod::MODE_CELL);
  mol.topo().update(AMod::CASE_CHANGE_MODE);
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 9;

  cout << "*********************************************" << endl;
  cout << "Check 7: change topo mode (cell)" << endl;
  cout << "*********************************************" << endl;
  mol.topo().mode(AMod::MODE_CELL);
  mol.topo().update(AMod::CASE_CHANGE_MODE);
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 10;

  cout << "*********************************************" << endl;
  cout << "Check 8: del # atoms" << endl;
  cout << "*********************************************" << endl;
  m = mol.natoms();
  n = m*(1.0-Util::RNG::uniform_std());
  cout << "# = " << n << endl;
  for(i = 0; i < n; i++) {
    mol.deleteAtom(mol.natoms()*Util::RNG::uniform_std());
  }
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 11;

  cout << "*********************************************" << endl;
  cout << "Check 9: new # atoms" << endl;
  cout << "*********************************************" << endl;
  n = m*(1.0-Util::RNG::uniform_std());
  cout << "# = " << n << endl;
  for(i = 0; i < n; i++) {
    data.type = AMod::PTE::Carbon;
    for(j = 0; j < 3; j++) data.x[j] = Util::RNG::uniform_std()*20.0;
    mol.newAtom(-1,data);
  }
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 12;

  cout << "*********************************************" << endl;
  cout << "Check 10: move # atoms" << endl;
  cout << "*********************************************" << endl;
  n = mol.natoms()*(1.0-Util::RNG::uniform_std());
  cout << "# = " << n << endl;
  for(i = 0; i < n; i++) {
    for(j = 0; j < 3; j++) data.x[j] = Util::RNG::uniform_std()*20.0;
    mol.moveAtom(mol.natoms()*Util::RNG::uniform_std(),data.x);
  }
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 13;

  cout << "*********************************************" << endl;
  cout << "Check 11: change axis" << endl;
  cout << "*********************************************" << endl;
  for(i = 0; i < mol.naxes(); i++)
    for(j = 0; j < 3; j++)
      mol.axis(i).x[j] += Util::RNG::uniform_std()*10.0;
  mol.topo().update(AMod::CASE_AXES);
  mol.topo().update(AMod::CASE_BRUTAL);
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 14;

  cout << "*********************************************" << endl;
  cout << "Check 12: adjust" << endl;
  cout << "*********************************************" << endl;
  madjuster.adjust();
  flag = mchecker.check();
  cout << endl;
  if(!flag) return 15;

  return 0;
}
