#include "MC_REBO.h"
#include "../Util/utilSTL.h"
#include <cmath>
#include <iostream>

namespace MC {

MC_REBO::MC_REBO() {}

MC_REBO::MC_REBO(AMod::Molecule& _mol): MCBasic(_mol) {}

MC_REBO::~MC_REBO() {}

void MC_REBO::init() {
  clear();
  AMod::Molecule& mol = molecule();
  mol.topo().bondRange() = REBO::Energy::bondRange();
  mol.topo().cellLength() = REBO::Energy::cellLength();
  computeEnergy();
}

void MC_REBO::final() { clear(); }

void MC_REBO::save() {
  molBak = molecule();
  energyBak = energy();
}

void MC_REBO::restore() {
  molecule() = molBak;
  energy() = energyBak;
}

double MC_REBO::computeEnergy() {
  AMod::Molecule& mol = molecule();
  mol.topo().update(AMod::CASE_AXES);
  mol.topo().update(AMod::CASE_GENERAL);
  mol.potential(energy().computeE(mol.topo()));
  return mol.potential();
}

double MC_REBO::dEnergyNewAtom(int _atomID, 
			       const AMod::Atom::Data& _atomData) {
  double dE;
  int k, kmax, kmax1, kmax2;
  AMod::Molecule& mol = molecule();
  AMod::Atom& a = mol.newAtom(_atomID, _atomData);
  //new bonds data
  for(k = 0, kmax = a.arrows.size(); k < kmax; k++)
    energy().bondData.push_back();
  //get 1st and 2nd neighbor
  bN1.clear(); bN2.clear();
  getNeighbors(_atomID);
  //old energy for affected bonds/atoms
  dE = 0.0;
  kmax1 = bN1.size();
  kmax2 = bN2.size();
  for(k = 0; k < kmax2; k++)
    dE -= energy().bondData[bN2[k]->id()].E;
  //bak
  potentialBak = mol.potential();
  atomIDBak = _atomID;
  bakAffected();
  //calculate bond energy for the neighboring bonds
  for(k = 0; k < kmax1; k++) energy().compute_r0(*bN1[k]);
  for(k = 0; k < kmax1; k++) energy().compute_r1(*bN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r1(*bN2[k]);
  for(k = 0; k < kmax1; k++) energy().compute_r2(*bN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r2(*bN2[k]);
  for(k = 0; k < kmax1; k++) dE += energy().bondData[bN1[k]->id()].E;  
  for(k = 0; k < kmax2; k++) dE += energy().bondData[bN2[k]->id()].E;
  mol.potential(mol.potential()+dE);
  return dE;
}

double MC_REBO::dEnergyDeleteAtom(int _atomID) {
  double dE;
  int k, kmax1, kmax2;
  AMod::Molecule& mol = molecule();
  AMod::Atom& a = mol.atom(_atomID);
  //get 1st and 2nd neighbor
  bN1.clear(); bN2.clear();
  getNeighbors(_atomID);
  //old energy for deleted/affected bonds/atoms
  dE = 0.0;
  kmax1 = bN1.size();
  kmax2 = bN2.size();
  for(k = 0; k < kmax1; k++)
    dE -= energy().bondData[bN1[k]->id()].E;
  for(k = 0; k < kmax2; k++)
    dE -= energy().bondData[bN2[k]->id()].E;
  //bak
  potentialBak = mol.potential();
  atomIDBak = _atomID;
  atomDataBak = a.data();
  bakDeleted();
  //delete atom & bonds
  for(k = 0; k < kmax1; k++) {
    energy().bondData.erase(bN1[k]->id());
    mol.topo().deleteBond(bN1[k]->id());
  }
  mol.deleteAtom(_atomID);
  //bak
  bakAffected();
  //calculate bond energy for the neighboring bonds
  for(k = 0; k < kmax2; k++) energy().compute_r1(*bN2[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r2(*bN2[k]);
  for(k = 0; k < kmax2; k++) dE += energy().bondData[bN2[k]->id()].E;
  mol.potential(mol.potential()+dE);
  return dE;
}

double MC_REBO::dEnergyMoveAtom(int _atomID, const double* _x) {
  double dE;
  int k, kmax, kmax1, kmax2;
  AMod::Molecule& mol = molecule();
  AMod::Atom& a = mol.atom(_atomID);
  //get 1st and 2nd neighbor
  bN1.clear(); bN2.clear();
  getNeighbors(_atomID);
  //old energy for deleted bonds/atoms
  dE = 0.0;
  kmax1 = bN1.size();
  for(k = 0; k < kmax1; k++)
    dE -= energy().bondData[bN1[k]->id()].E;
  //bak
  potentialBak = mol.potential();
  atomIDBak = _atomID;
  for(k = 0; k < 3; k++) xBak[k] = a->x[k];
  bakDeleted();
  //delete bonds data
  for(k = 0; k < kmax1; k++) {
    energy().bondData.erase(bN1[k]->id());
    mol.topo().deleteBond(bN1[k]->id());
  }
  //move atom
  mol.moveAtom(_atomID,_x);
  //new bonds data
  for(k = 0, kmax = a.arrows.size(); k < kmax; k++)
    energy().bondData.push_back();
  //get 1st and 2nd neighbor
  bN1.clear();
  getNeighbors(_atomID);
  //old energy for affected bonds/atoms
  kmax1 = bN1.size();
  kmax2 = bN2.size();
  for(k = 0; k < kmax2; k++)
    dE -= energy().bondData[bN2[k]->id()].E;
  //bak
  bakAffected();
  //calculate bond energy for the neighboring bonds
  for(k = 0; k < kmax1; k++) energy().compute_r0(*bN1[k]);
  for(k = 0; k < kmax1; k++) energy().compute_r1(*bN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r1(*bN2[k]);
  for(k = 0; k < kmax1; k++) energy().compute_r2(*bN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r2(*bN2[k]);
  for(k = 0; k < kmax1; k++) dE += energy().bondData[bN1[k]->id()].E;  
  for(k = 0; k < kmax2; k++) dE += energy().bondData[bN2[k]->id()].E;
  mol.potential(mol.potential()+dE);
  return dE;
}

void MC_REBO::unNewAtom() {
  restoreAffected();
  int k, kmax;
  AMod::Atom& a = molecule().atom(atomIDBak);
  for(k = 0, kmax = a.arrows.size(); k < kmax; k++)
    energy().bondData.pop_back();
  molecule().deleteAtom(atomIDBak);
  molecule().potential(potentialBak);
}

void MC_REBO::unDeleteAtom() {
  restoreAffected();
  molecule().newLonelyAtom(atomIDBak,atomDataBak);
  restoreDeleted();
  molecule().potential(potentialBak);
}

void MC_REBO::unMoveAtom() {
  restoreAffected();
  int k, kmax;
  AMod::Atom& a = molecule().atom(atomIDBak);
  for(k = 0, kmax = a.arrows.size(); k < kmax; k++)
    energy().bondData.pop_back();
  molecule().makeLonelyAtom(atomIDBak);
  molecule().moveLonelyAtom(atomIDBak,xBak);
  restoreDeleted();
  molecule().potential(potentialBak);
}

void MC_REBO::clear() {
  bN1.clear();
  bN2.clear();
  molBak.reset();
  energyBak.reset();
  atomDataBak.reset();
  bondDataIDBak.clear();
  bondDataBak.clear();
  bondBak.clear();
  bondDataDelBak.clear();
}

void MC_REBO::bakAffected() {
  int k, kmax = bN2.size();
  bondDataIDBak.resize(kmax);
  bondDataBak.resize(kmax);
  for(k = 0; k < kmax; k++) {
    bondDataIDBak[k] = bN2[k]->id();
    bondDataBak[k] = energy().bondData[bN2[k]->id()];
  }
}

void MC_REBO::bakDeleted() {
  int k, kmax = bN1.size();
  const int* n;
  bondBak.resize(5*kmax);
  bondDataDelBak.resize(kmax);
  for(k = 0; k < kmax; k++) {
    n = bN1[k]->arrow(0).n;
    bondBak[5*k  ] = bN1[k]->arrow(0).host().id();
    bondBak[5*k+1] = bN1[k]->arrow(1).host().id();
    bondBak[5*k+2] = n[0];
    bondBak[5*k+3] = n[1];
    bondBak[5*k+4] = n[2];
    bondDataDelBak[k] = energy().bondData[bN1[k]->id()];
  }
}

void MC_REBO::restoreAffected() {
  int k, kmax = bondDataIDBak.size();
  for(k = 0; k < kmax; k++) {
    energy().bondData[bondDataIDBak[k]] = bondDataBak[k];
  }
}

void MC_REBO::restoreDeleted() {
  AMod::Molecule& mol = molecule();
  int n[3], k, kmax = bondBak.size()/5;
  for(k = 0; k < kmax; k++) {
    n[0] = bondBak[5*k+2];
    n[1] = bondBak[5*k+3];
    n[2] = bondBak[5*k+4];
    mol.topo().newBond(bondBak[5*k], bondBak[5*k+1], n);
    energy().bondData.push_back(bondDataDelBak[k]);
  }
}

void MC_REBO::getNeighbors(int _atomID) {
  int i, imax, j, jmax, k, kmax;
  AMod::Molecule& mol = molecule();
  AMod::Atom& a = mol.atom(_atomID);
  //bN1
  for(k = 0, kmax = a.arrows.size(); k < kmax; k++) {
    Util::push_back_unique(bN1, &a.arrows[k].bond());
  }
  //bN2
  for(k = 0, kmax = a.arrows.size(); k < kmax; k++) {
    AMod::Atom& a2 = a.arrows[k].moon();
    for(j = 0, jmax = a2.arrows.size(); j < jmax; j++) {
      AMod::Atom& a3 = a2.arrows[j].moon();
      if(a3.id() == a.id()) continue;
      for(i = 0, imax = a3.arrows.size(); i < imax; i++) {
	if(a3.arrows[i].moon().id() != a.id())
	  Util::push_back_unique(bN2, &a3.arrows[i].bond());
      }//end for(i...
    }//end for(j...
  }//end for(k...
}

AMod::MolTopo& MC_REBO::molTopo() { return molecule().topo(); }

const AMod::MolTopo& MC_REBO::molTopo() const { return molecule().topo(); }

}/* MC */
