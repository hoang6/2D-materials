#include "MolChecker.h"
#include "Molecule.h"
#include <iostream>
#include <algorithm>
#include <set>

namespace AMod {

MolChecker::MolChecker() {}

MolChecker::MolChecker(const Molecule& _mol): ConstMolServ(_mol) {}

void MolChecker::reset() {}

bool MolChecker::neighbCheck() const{
  using namespace std;
  cout << "[Check Neighb" << endl;
  set<int> pool;
  const Molecule& mol = molecule();
  for(int k = 0, kmax = mol.natoms(); k < kmax; k++) {
    pool.clear();
    const Atom& a = mol.atom(k);
    for(int j = 0, jmax = a.arrows.size(); j < jmax; j++) {
      if(!pool.insert(a.arrows[j].moon().id()).second) {
	cout << "ERROR: atom " << k << " has duplicate neighbors "
	     << a.arrows[j].moon().id() << endl;
	return false;
      }
    }//end for(int j...
  }//end for(int k...
  cout << "]" << endl;
  return true;
}

bool MolChecker::connectionCheck() const {
  using namespace std;
  cout << "[Check Connection" << endl;
  const Molecule& mol = molecule();
  int natoms = mol.natoms();
  int nbonds = mol.nbonds();
  int npairs = (natoms*(natoms-1))/2;
  std::vector<bool> visited(npairs, false);
  int i, j, k, index;
  for(k = 0; k < nbonds; k++) {
    const Bond& bond = mol.bond(k);
    const Atom& a0 = bond.arrow(0).host();
    const Atom& a1 = bond.arrow(1).host();
    i = a0.id();
    j = a1.id();
    if(i == j) return false;
    if(i > j) std::swap(i, j);
    index = natoms*i-((i+1)*(i+2))/2+j;
    if(visited[index]) continue;
    if(!mol.topo().bondok(k)) {
      cout << "ERROR: bond " << k << " is not in bondRange" << endl; 
      //debug
      printAtom(i);
      printAtom(j);
      return false;
    }
    visited[index] = true;
  }//end for(k = 0; k < nbonds; k++)
  for(i = 0; i < natoms; i++) 
    for(j = i+1; j < natoms; j++) {
      index = natoms*i-((i+1)*(i+2))/2+j;
      if(visited[index]) continue;
      if(mol.topo().bondRange().ok(mol.atom(i)->type,
				   mol.atom(j)->type,
				   mol.topo().distance(i,j))) {
	cout << "ERROR: "
	     << "atom " << i << " and " 
	     << "atom " << j << " should be bonded" << endl;
	//debug
	cout << "dist = " << mol.topo().distance(i,j) << endl;
	printAtom(i);
	printAtom(j);
	return false;
      }
      visited[index] = true;
    }
  cout << "]" << endl;
  return true;
}

bool MolChecker::topoCheck() const {
  using namespace std;
  cout << "[Check Topo" << endl;
  const MolTopo& topo = molecule().topo();  
  if(topo.mode() == MODE_NULL) {
    cout << "mode = MODE_NONE" << endl; 
  }
  else if(topo.mode() == MODE_FREEZE) {
    cout << "mode = MODE_FREEZE" << endl; 
  }
  else if(topo.mode() == MODE_CELL) {
    cout << "mode = MODE_CELL" << endl;    
  }
  cout << "<Print Cell>" << endl;
  cout << "cellLength = " << topo.cellLength() << endl;
  cout << "cellAxes = ";
  for(int i = 0; i < 3; i++) {
    printVec(topo.cellAxes[i],3); cout << " ";
  }
  cout << endl;
  cout << "ncell = "; printVec(topo.ncell,3); cout << endl;
  cout << "natoms = " << topo.natoms() << endl;
  cout << "nbonds = " << topo.nbonds() << endl;
  cout << "c2aLinks.size = " << topo.c2aLinks.size() << endl;
  cout << "a2cLinks.size = " << topo.a2cLinks.size() << endl;
  if(topo.mode() == MODE_CELL) {
    MolTopo::Cell tc;
    /* check c2aLinks */
    cout << "<Check c2aLinks" << endl;
    MolTopo::C2ALinks::const_iterator cit;
    MolTopo::AtomPtrList::const_iterator ait;
    bool flag;
    for(int k = 0, kmax = topo.natoms(); k < kmax; k++) {
      topo.whatCell(topo.atom(k)->x,tc);
      cit = topo.c2aLinks.find(tc);
      flag = false;
      if(cit != topo.c2aLinks.end()) {
	const MolTopo::AtomPtrList& apl = (*cit).second;
	for(ait = apl.begin(); ait != apl.end(); ait++) {
	  if((*ait)->id() == k) { flag = true; break; }
	}
      }
      if(!flag) {
	cout << "ERROR: atom " << k << " not in c2aLinks" << endl;
	return false;
      }
    }//end for(int k...
    cout << ">" << endl;
    /* check a2cLinks */
    cout << "<Check a2cLinks" << endl;
    if(topo.natoms() != topo.a2cLinks.size()) {
      cout << "ERROR: natoms != a2cLinks.size" << endl;
      return false;
    }
    for(int k = 0, kmax = topo.natoms(); k < kmax; k++) {
      topo.whatCell(topo.atom(k)->x,tc);
      if(!(topo.a2cLinks[k].cell == tc)) {
	cout << "ERROR: atom " << k << " has wrong cell, i.e. ";
	printVec(topo.a2cLinks[k].cell.x,3); cout << " != "; printVec(tc.x,3);
	cout << endl;
	return false;
      }
    }//end for(k...
    cout << ">" << endl;
  }//end if(topo().mode()...
  cout << "]" << endl;
  return true;
}

bool MolChecker::check() const {
  int k;
  bool flag;
  bool flags[] = {
    neighbCheck(), 
    connectionCheck(), 
    topoCheck()
  };
  std::cout << "---------------------------------------------" << std::endl;
  showMessage("neighbCheck", flags[0]);
  showMessage("connectionCheck", flags[1]);
  showMessage("topoCheck", flags[2]);
  std::cout << "---------------------------------------------" << std::endl;
  flag = true;
  for(k = 0; k < int(sizeof(flags)/sizeof(bool)); k++) flag = (flag&&flags[k]);
  return flag;
}

void MolChecker::printAtom(int _atomID) const {
  using namespace std;
  const Molecule& mol = molecule();
  const Atom& a = mol.atom(_atomID);
  cout << "@atom " << a.id() << endl;
  cout << "type = " << a->type << endl;
  cout << "x = "; printVec(a->x,3); cout << endl;
  if(mol.topo().mode() == MODE_CELL) {
    MolTopo::Cell tc;
    mol.topo().whatCell(a->x,tc);
    cout << "cell = "; printVec(tc.x,3); cout << endl;
  }
  cout << "@" << endl;
}

template<class T>
void MolChecker::printVec(const T* v, long len) {
  std::cout << "{";
  for(long i = 0; i < len; i++) 
    std::cout << v[i] << (i==len-1 ? "" : ",");
  std::cout << "}";
}

void MolChecker::showMessage(const char* whatCheck, bool passed) const {
  std::cout << whatCheck << " " << (passed ? "passed" : "failed") << std::endl;
}

} /* AMod */
