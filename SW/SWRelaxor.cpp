#include "SWRelaxor.h"
#include "../AMod/Molecule.h"
#include "../AMod/PTE.h"
#include "../AMod/constants.h"
#include "../Util/FuncDFunc.h"
#include "../Util/ConjGrad.h"
#include <set>
#include <cmath>

namespace SW {

SWRelaxor::SWRelaxor() {}

SWRelaxor::SWRelaxor(AMod::Molecule& _mol): AMod::MolServ(_mol) {} 

double SWRelaxor::relax(int _bondID, double tol) {
  setup(_bondID);
  int k, iter, n = nDOF();
  double fmin;
  AMod::Molecule& mol = molecule();
  double* _x = new double[n];

  for(k = 0; k < n; k++) 
    _x[k] = mol.atom(atomIDs[k/3])->x[k%3];
  Util::FuncDFunc<SWRelaxor> fdf(*this, n, 
				 &SWRelaxor::func, 
				 &SWRelaxor::dfunc);
  Util::conjGrad(_x, n, tol, 10000, fdf, iter, fmin);

  delete [] _x;
  return fmin;
}

void SWRelaxor::setup(int _bondID) {
  int j, k, jmax, kmax;
  std::set<int> ids;
  const AMod::Molecule& mol = molecule();
  atomIDs.clear();
  ids.clear();
  for(k = 0; k < 2; k++) {
    const AMod::Atom& a = mol.bond(_bondID).arrow(k).host();
    jmax = a.arrows.size();
    for(j = 0; j < jmax; j++) {
      const AMod::Arrow& ar = a.arrows[j];
      if(ids.insert(ar.moon().id()).second)
	atomIDs.push_back(ar.moon().id());
    }//end for(j...
  }//end for(k...
  atomIDPairs.clear();
  ids.clear();
  kmax = atomIDs.size();
  for(k = 0; k < kmax; k++) {
    const AMod::Atom& a = mol.atom(atomIDs[k]);
    jmax = a.arrows.size();
    for(j = 0; j < jmax; j++) {
      if(ids.insert(a.arrows[j].bond().id()).second)
	atomIDPairs.push_back(AMod::DInt(a.id(),a.arrows[j].moon().id()));
    }
  }
}

int SWRelaxor::nDOF() const { return 3*atomIDs.size(); }

double SWRelaxor::func(double* _x) {
  int k, id1, id2;
  int natoms = atomIDs.size();
  int npairs = atomIDPairs.size();
  AMod::Molecule& mol = molecule();
  double rref, dr, score;
  for(k = 0; k < natoms; k++)
    mol.moveAtom(atomIDs[k], _x+3*k);
  score = 0.0;
  for(k = 0; k < npairs; k++) {
    id1 = atomIDPairs[k].first;
    id2 = atomIDPairs[k].second;
    const AMod::Atom& a1 = mol.atom(id1);
    const AMod::Atom& a2 = mol.atom(id2);
    const AMod::DInt idpair = atomIDPairs[k];
    rref = stdBondLength(a1->type,a2->type);
    dr = mol.topo().distance(id1,id2)-rref;
    score += dr*dr;
  }
  return score;
}

void SWRelaxor::dfunc(double* _x, double* _g) {
  int k, n = nDOF();
  double score1, score2, temp;
  const double step = 1e-6;
  for(k = 0; k < n; k++) {
    temp = _x[k];
    _x[k] += step;
    score1 = func(_x);
    _x[k] -= 2.0*step;
    score2 = func(_x);
    _g[k] = (score1-score2)/2.0/step;
    _x[k] = temp;
  }
}

double SWRelaxor::stdBondLength(int type1, int type2) {
  if(type1 == AMod::PTE::Carbon && type2 == AMod::PTE::Carbon) 
    return AMod::CC_BOND_LENGTH;
  if(type1 == AMod::PTE::Carbon && type2 == AMod::PTE::Hydrogen) 
    return AMod::CH_BOND_LENGTH;
  if(type1 == AMod::PTE::Hydrogen && type2 == AMod::PTE::Carbon) 
    return AMod::CH_BOND_LENGTH;
  if(type1 == AMod::PTE::Hydrogen && type2 == AMod::PTE::Hydrogen) 
    return AMod::HH_BOND_LENGTH;
  return 0.0;
}

double SWRelaxor::stdBondLength(const AMod::Bond& bond) {
  return stdBondLength(bond.arrow(0).host()->type, bond.arrow(1).host()->type);
}

}/* SW */
