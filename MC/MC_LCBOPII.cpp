#include "MC_LCBOPII.h"
#include "../Util/utilSTL.h"
#include <cmath>
#include <iostream>

namespace MC {

MC_LCBOPII::MC_LCBOPII() {}

MC_LCBOPII::MC_LCBOPII(AMod::Molecule& _mol): MCBasic(_mol) {}

MC_LCBOPII::~MC_LCBOPII() {}

void MC_LCBOPII::init() {
  clear();
  AMod::Molecule& mol = molecule();
  //SR
  mol.topo().bondRange() = LCBOPII::Energy::bondRangeSR();
  mol.topo().cellLength() = LCBOPII::Energy::cellLengthSR();
  mol.topo().update(AMod::CASE_CHANGE_MODE);
  //MR
  topoMR.reset();
  topoMR.attach(molecule());
  topoMR.sync();
  topoMR.mode(mol.topo().mode());
  topoMR.bondRange() = LCBOPII::Energy::bondRangeMR();
  topoMR.cellLength() = LCBOPII::Energy::cellLengthMR();
  topoMR.update(AMod::CASE_CHANGE_MODE);
  //LR
  topoLR.reset();
  topoLR.attach(molecule());
  topoLR.sync();
  topoLR.mode(mol.topo().mode());
  topoLR.bondRange() = LCBOPII::Energy::bondRangeLR();
  topoLR.cellLength() = LCBOPII::Energy::cellLengthLR();
  topoLR.update(AMod::CASE_CHANGE_MODE);
  //energy
  computeEnergy();
}

void MC_LCBOPII::final() { clear(); }

void MC_LCBOPII::save() {
  molBak = molecule();
  topoMRBak = topoMR;
  topoLRBak = topoLR;
  energyBak = energy();
}

void MC_LCBOPII::restore() {
  molecule() = molBak;
  topoMR = topoMRBak;
  topoLR = topoLRBak;
  energy() = energyBak;
}

double MC_LCBOPII::computeEnergy() {
  AMod::Molecule& mol = molecule();
  mol.topo().update(AMod::CASE_AXES);
  mol.topo().update(AMod::CASE_GENERAL);
  topoMR.update(AMod::CASE_AXES);
  topoMR.update(AMod::CASE_GENERAL);
  topoLR.update(AMod::CASE_AXES);
  topoLR.update(AMod::CASE_GENERAL);
  mol.potential(energy().computeE(mol.topo(),topoMR,topoLR));
  return mol.potential();
}

double MC_LCBOPII::dEnergyNewAtom(int _atomID, 
				  const AMod::Atom::Data& _atomData) {
  double dE;
  int k, kmax, kmax0, kmax1, kmax2, kmax3, kmax4;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.newAtom(_atomID, _atomData);
  AMod::Atom& aLR = topoLR.newAtom(_atomID, aSR.data());
  AMod::Atom& aMR = topoMR.newAtom(_atomID, aSR.data());
  //new bonds/atoms data
  for(k = 0, kmax = aSR.arrows.size(); k < kmax; k++)
    energy().bondDataSR.push_back();
  for(k = 0, kmax = aLR.arrows.size(); k < kmax; k++)
    energy().bondDataLR.push_back();
  energy().atomDataMR.insert(_atomID);
  //get 1st and 2nd neighbors
  bSRN1.clear(); bSRN2.clear(); 
  aSRN1.clear(); aMRN2.clear();
  getNeighbors(_atomID);
  //old energy for affected bonds/atoms
  dE = 0.0;
  kmax0 = aSRN1.size();
  kmax1 = bSRN1.size();
  kmax2 = bSRN2.size();
  kmax3 = aLR.arrows.size();
  kmax4 = aMRN2.size();
  for(k = 0; k < kmax2; k++) 
    dE -= energy().bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax4; k++)
    dE -= energy().atomDataMR[aMRN2[k]->id()].E;
  //bak
  potentialBak = mol.potential();
  atomIDBak = _atomID;
  bakAffected();
  //new energy for new/affected bonds/atoms
  for(k = 0; k < kmax1; k++) energy().compute_r0(*bSRN1[k]);
  for(k = 0; k < kmax0; k++) energy().compute_r1(*aSRN1[k]);
  energy().compute_r1(aSR);
  for(k = 0; k < kmax1; k++) energy().compute_r2(*bSRN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r2(*bSRN2[k]);
  for(k = 0; k < kmax3; k++) energy().compute_r3(aLR.arrows[k].bond());
  for(k = 0; k < kmax4; k++) energy().compute_r4(*aMRN2[k],mol.topo());
  energy().compute_r4(aMR,mol.topo());
  //sum bonds/atoms energy
  for(k = 0; k < kmax1; k++) 
    dE += energy().bondDataSR[bSRN1[k]->id()].E;
  for(k = 0; k < kmax2; k++) 
    dE += energy().bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax3; k++) 
    dE += energy().bondDataLR[aLR.arrows[k].bond().id()].E;
  for(k = 0; k < kmax4; k++)
    dE += energy().atomDataMR[aMRN2[k]->id()].E;
  dE += energy().atomDataMR[aMR.id()].E;
  mol.potential(mol.potential()+dE);
  return dE;  
}

double MC_LCBOPII::dEnergyDeleteAtom(int _atomID) {
  double dE;
  int k, kmax0, kmax1, kmax2, kmax3, kmax4;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(_atomID);
  AMod::Atom& aLR = topoLR.atom(_atomID);
  AMod::Atom& aMR = topoMR.atom(_atomID);
  //get 1st and 2nd neighbors
  bSRN1.clear(); bSRN2.clear(); 
  aSRN1.clear(); aMRN2.clear();
  getNeighbors(_atomID);
  //old energy for deleted/affected bonds/atoms
  dE = 0.0;
  kmax0 = aSRN1.size();
  kmax1 = bSRN1.size();
  kmax2 = bSRN2.size();
  kmax3 = aLR.arrows.size();
  kmax4 = aMRN2.size();
  for(k = 0; k < kmax1; k++) 
    dE -= energy().bondDataSR[bSRN1[k]->id()].E;
  for(k = 0; k < kmax2; k++) 
    dE -= energy().bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax3; k++) 
    dE -= energy().bondDataLR[aLR.arrows[k].bond().id()].E;
  for(k = 0; k < kmax4; k++)
    dE -= energy().atomDataMR[aMRN2[k]->id()].E;
  dE -= energy().atomDataMR[aMR.id()].E;
  //bak
  potentialBak = mol.potential();
  atomIDBak = _atomID;
  atomDataBak = aSR.data();
  bakDeleted();
  //delete bonds data
  for(k = 0; k < kmax1; k++) {
    energy().bondDataSR.erase(bSRN1[k]->id());
    mol.topo().deleteBond(bSRN1[k]->id());
  }
  for(k = 0; k < kmax3; k++) {
    energy().bondDataLR.erase(aLR.arrows[0].bond().id());
    topoLR.deleteBond(aLR.arrows[0].bond().id());
  }
  energy().atomDataMR.erase(aMR.id());
  mol.deleteAtom(_atomID);
  topoLR.deleteAtom(_atomID);  
  topoMR.deleteAtom(_atomID);
  //bak
  bakAffected();
  //update energy for affected bonds/atoms
  for(k = 0; k < kmax0; k++) energy().compute_r1(*aSRN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r2(*bSRN2[k]);
  for(k = 0; k < kmax4; k++) energy().compute_r4(*aMRN2[k],mol.topo());
  //sum bonds/atoms energy
  for(k = 0; k < kmax2; k++) 
    dE += energy().bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax4; k++)
    dE += energy().atomDataMR[aMRN2[k]->id()].E;
  mol.potential(mol.potential()+dE);
  return dE;
}

double MC_LCBOPII::dEnergyMoveAtom(int _atomID, const double* _x) {
  double dE;
  int k, kmax, kmax0, kmax1, kmax2, kmax3, kmax4;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(_atomID);
  AMod::Atom& aLR = topoLR.atom(_atomID);
  AMod::Atom& aMR = topoMR.atom(_atomID);
  //get 1st and 2nd neighbors
  bSRN1.clear(); bSRN2.clear(); 
  aSRN1.clear(); aMRN2.clear();
  getNeighbors(_atomID);
  //old energy for deleted bonds/atoms
  dE = 0.0;
  kmax1 = bSRN1.size();
  kmax3 = aLR.arrows.size();
  for(k = 0; k < kmax1; k++) 
    dE -= energy().bondDataSR[bSRN1[k]->id()].E;
  for(k = 0; k < kmax3; k++) 
    dE -= energy().bondDataLR[aLR.arrows[k].bond().id()].E;
  dE -= energy().atomDataMR[aMR.id()].E;
  //bak
  potentialBak = mol.potential();
  atomIDBak = _atomID;
  for(k = 0; k < 3; k++) xBak[k] = aSR->x[k];
  bakDeleted();
  //delete bonds data
  for(k = 0; k < kmax1; k++) {
    energy().bondDataSR.erase(bSRN1[k]->id());
    mol.topo().deleteBond(bSRN1[k]->id());
  }
  for(k = 0; k < kmax3; k++) {
    energy().bondDataLR.erase(aLR.arrows[0].bond().id());
    topoLR.deleteBond(aLR.arrows[0].bond().id());
  }
  //move atom
  mol.moveAtom(_atomID,_x);
  topoMR.moveAtom(_atomID);
  topoLR.moveAtom(_atomID);
  //new bonds data
  for(k = 0, kmax = aSR.arrows.size(); k < kmax; k++)
    energy().bondDataSR.push_back();
  for(k = 0, kmax = aLR.arrows.size(); k < kmax; k++)
    energy().bondDataLR.push_back();
  //get 1st and 2nd neighbors
  bSRN1.clear(); 
  getNeighbors(_atomID);
  //old energy for affected bonds/atoms
  kmax0 = aSRN1.size();
  kmax1 = bSRN1.size();
  kmax2 = bSRN2.size();
  kmax3 = aLR.arrows.size();
  kmax4 = aMRN2.size();
  for(k = 0; k < kmax2; k++) 
    dE -= energy().bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax4; k++)
    dE -= energy().atomDataMR[aMRN2[k]->id()].E;
  //bak
  bakAffected();
  //new energy for new/affected bonds/atoms
  for(k = 0; k < kmax1; k++) energy().compute_r0(*bSRN1[k]);
  for(k = 0; k < kmax0; k++) energy().compute_r1(*aSRN1[k]);
  energy().compute_r1(aSR);
  for(k = 0; k < kmax1; k++) energy().compute_r2(*bSRN1[k]);
  for(k = 0; k < kmax2; k++) energy().compute_r2(*bSRN2[k]);
  for(k = 0; k < kmax3; k++) energy().compute_r3(aLR.arrows[k].bond());
  for(k = 0; k < kmax4; k++) energy().compute_r4(*aMRN2[k],mol.topo());
  energy().compute_r4(aMR,mol.topo());
  //sum bonds/atoms energy
  for(k = 0; k < kmax1; k++) 
    dE += energy().bondDataSR[bSRN1[k]->id()].E;
  for(k = 0; k < kmax2; k++) 
    dE += energy().bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax3; k++) 
    dE += energy().bondDataLR[aLR.arrows[k].bond().id()].E;
  for(k = 0; k < kmax4; k++)
    dE += energy().atomDataMR[aMRN2[k]->id()].E;
  dE += energy().atomDataMR[aMR.id()].E;
  mol.potential(mol.potential()+dE);
  return dE;
}

void MC_LCBOPII::unNewAtom() {
  restoreAffected();
  int k, kmax;
  AMod::Atom& aSR = molecule().atom(atomIDBak);
  AMod::Atom& aLR = topoLR.atom(atomIDBak);
  for(k = 0, kmax = aSR.arrows.size(); k < kmax; k++)
    energy().bondDataSR.pop_back();
  for(k = 0, kmax = aLR.arrows.size(); k < kmax; k++)
    energy().bondDataLR.pop_back();
  energy().atomDataMR.erase(atomIDBak);
  molecule().deleteAtom(atomIDBak);
  topoMR.deleteAtom(atomIDBak);
  topoLR.deleteAtom(atomIDBak);
  molecule().potential(potentialBak);
}

void MC_LCBOPII::unDeleteAtom() {
  restoreAffected();
  AMod::Atom& aSR = molecule().newLonelyAtom(atomIDBak,atomDataBak);
  topoLR.newLonelyAtom(atomIDBak,aSR.data());
  topoMR.newLonelyAtom(atomIDBak,aSR.data());
  energy().atomDataMR.insert(atomIDBak);
  restoreDeleted();
  molecule().potential(potentialBak);
}

void MC_LCBOPII::unMoveAtom() {
  restoreAffected();
  int k, kmax;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(atomIDBak);
  AMod::Atom& aLR = topoLR.atom(atomIDBak);
  for(k = 0, kmax = aSR.arrows.size(); k < kmax; k++)
    energy().bondDataSR.pop_back();
  for(k = 0, kmax = aLR.arrows.size(); k < kmax; k++)
    energy().bondDataLR.pop_back();
  mol.topo().makeLonelyAtom(atomIDBak);
  topoMR.makeLonelyAtom(atomIDBak);
  topoLR.makeLonelyAtom(atomIDBak);
  mol.moveLonelyAtom(atomIDBak,xBak);
  topoMR.moveLonelyAtom(atomIDBak);
  topoLR.moveLonelyAtom(atomIDBak);
  restoreDeleted();
  molecule().potential(potentialBak);
}

void MC_LCBOPII::clear() {
  bSRN1.clear();
  bSRN2.clear();
  aSRN1.clear();
  aMRN2.clear();
  molBak.reset();
  topoMRBak.reset();
  topoLRBak.reset();
  energyBak.reset();
  atomDataBak.reset();
  bondDataSRIDBak.clear();
  atomDataMRIDBak.clear();
  bondDataSRBak.clear();
  atomDataMRBak.clear();
  bondSRBak.clear();
  bondLRBak.clear();
  bondMRBak.clear();
  bondDataSRDelBak.clear();
  bondDataLRDelBak.clear();
}

void MC_LCBOPII::bakAffected() {
  int k;
  int kmax2 = bSRN2.size();
  int kmax4 = aMRN2.size();
  bondDataSRIDBak.resize(kmax2);
  bondDataSRBak.resize(kmax2);
  for(k = 0; k < kmax2; k++) {
    bondDataSRIDBak[k] = bSRN2[k]->id();
    bondDataSRBak[k] = energy().bondDataSR[bSRN2[k]->id()];
  }
  atomDataMRIDBak.resize(kmax4);
  atomDataMRBak.resize(kmax4);
  for(k = 0; k < kmax4; k++) {
    atomDataMRIDBak[k] = aMRN2[k]->id();
    atomDataMRBak[k] = energy().atomDataMR[aMRN2[k]->id()];
  }
}

void MC_LCBOPII::bakDeleted() {
  int k;
  const int* n;
  AMod::Atom& aLR = topoLR.atom(atomIDBak);
  AMod::Atom& aMR = topoMR.atom(atomIDBak);
  int kmax1 = bSRN1.size();
  int kmax3 = aLR.arrows.size();
  int kmax5 = aMR.arrows.size();
  bondSRBak.resize(5*kmax1);
  bondDataSRDelBak.resize(kmax1);
  for(k = 0; k < kmax1; k++) {
    n = bSRN1[k]->arrow(0).n;
    bondSRBak[5*k  ] = bSRN1[k]->arrow(0).host().id();
    bondSRBak[5*k+1] = bSRN1[k]->arrow(1).host().id();
    bondSRBak[5*k+2] = n[0];
    bondSRBak[5*k+3] = n[1];
    bondSRBak[5*k+4] = n[2];
    bondDataSRDelBak[k] = energy().bondDataSR[bSRN1[k]->id()];
  }
  bondLRBak.resize(5*kmax3);
  bondDataLRDelBak.resize(kmax3);
  for(k = 0; k < kmax3; k++) {
    n = aLR.arrows[k].n;
    bondLRBak[5*k  ] = aLR.id();
    bondLRBak[5*k+1] = aLR.arrows[k].moon().id();
    bondLRBak[5*k+2] = n[0];
    bondLRBak[5*k+3] = n[1];
    bondLRBak[5*k+4] = n[2];
    bondDataLRDelBak[k] = energy().bondDataLR[aLR.arrows[k].bond().id()];
  }
  bondMRBak.resize(5*kmax5);
  for(k = 0; k < kmax5; k++) {
    n = aMR.arrows[k].n;
    bondMRBak[5*k  ] = aMR.id();
    bondMRBak[5*k+1] = aMR.arrows[k].moon().id();
    bondMRBak[5*k+2] = n[0];
    bondMRBak[5*k+3] = n[1];
    bondMRBak[5*k+4] = n[2];
  }
  atomDatumMRDelBak = energy().atomDataMR[aMR.id()];
}

void MC_LCBOPII::restoreAffected() {
  int k, kmax;
  for(k = 0, kmax = bondDataSRIDBak.size(); k < kmax; k++) {
    energy().bondDataSR[bondDataSRIDBak[k]] = bondDataSRBak[k];
  }
  for(k = 0, kmax = atomDataMRIDBak.size(); k < kmax; k++) {
    energy().atomDataMR[atomDataMRIDBak[k]] = atomDataMRBak[k];
  }
}

void MC_LCBOPII::restoreDeleted() {
  AMod::Molecule& mol = molecule();
  int n[3], k, kmax;
  for(k = 0, kmax = bondSRBak.size()/5; k < kmax; k++) {
    n[0] = bondSRBak[5*k+2];
    n[1] = bondSRBak[5*k+3];
    n[2] = bondSRBak[5*k+4];
    mol.topo().newBond(bondSRBak[5*k], bondSRBak[5*k+1], n);
    energy().bondDataSR.push_back(bondDataSRDelBak[k]);
  }
  for(k = 0, kmax = bondLRBak.size()/5; k < kmax; k++) {
    n[0] = bondLRBak[5*k+2];
    n[1] = bondLRBak[5*k+3];
    n[2] = bondLRBak[5*k+4];
    topoLR.newBond(bondLRBak[5*k], bondLRBak[5*k+1], n);
    energy().bondDataLR.push_back(bondDataLRDelBak[k]);
  }
  for(k = 0, kmax = bondMRBak.size()/5; k < kmax; k++) {
    n[0] = bondMRBak[5*k+2];
    n[1] = bondMRBak[5*k+3];
    n[2] = bondMRBak[5*k+4];
    topoMR.newBond(bondMRBak[5*k], bondMRBak[5*k+1], n);
  }
  energy().atomDataMR[atomIDBak] = atomDatumMRDelBak;
}

void MC_LCBOPII::getNeighbors(int _atomID) {
  int i, imax, j, jmax, k, kmax;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(_atomID);
  AMod::Atom& aMR = topoMR.atom(_atomID);
  //bSRN1
  for(k = 0, kmax = aSR.arrows.size(); k < kmax; k++) {
    Util::push_back_unique(bSRN1, &aSR.arrows[k].bond());
  }
  //aMRN2
  for(k = 0, kmax = aMR.arrows.size(); k < kmax; k++) {
    Util::push_back_unique(aMRN2, &aMR.arrows[k].moon());
  }
  //affected bonds/atoms (bSRN2 & aMRN2)
  for(k = 0, kmax = aSR.arrows.size(); k < kmax; k++) {
    AMod::Atom& aSR2 = aSR.arrows[k].moon();
    //aSRN1
    Util::push_back_unique(aSRN1, &aSR2);
    //aMRN2
    Util::push_back_unique(aMRN2, &topoMR.atom(aSR2.id()));
    for(j = 0, jmax = aSR2.arrows.size(); j < jmax; j++) {
      AMod::Atom& aSR3 = aSR2.arrows[j].moon();
      if(aSR3.id() == aSR.id()) continue;
      //aMRN2
      Util::push_back_unique(aMRN2, &topoMR.atom(aSR3.id()));
      for(i = 0, imax = aSR3.arrows.size(); i < imax; i++) {
	AMod::Atom& aSR4 = aSR3.arrows[i].moon();
	if(aSR4.id() != aSR.id()) {
	  //bSRN2
	  Util::push_back_unique(bSRN2, &aSR3.arrows[i].bond());
	  //aMRN2
	  if(aSR4.id() != aSR2.id())
	    Util::push_back_unique(aMRN2, &topoMR.atom(aSR4.id()));
	}
      }//end for(i = 0...
    }//end for(j = 0...
  }//end for(k  = 0...
}

AMod::MolTopo& MC_LCBOPII::molTopoSR() { return molecule().topo(); }

AMod::MolTopo& MC_LCBOPII::molTopoMR() { return topoMR; }

AMod::MolTopo& MC_LCBOPII::molTopoLR() { return topoLR; }

const AMod::MolTopo& MC_LCBOPII::molTopoSR() const { return molecule().topo(); }

const AMod::MolTopo& MC_LCBOPII::molTopoMR() const { return topoMR; }

const AMod::MolTopo& MC_LCBOPII::molTopoLR() const { return topoLR; }

}/* MC */
