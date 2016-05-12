#include "KernelLCBOPIIN.h"
#include "Molecule.h"

namespace AMod {

KernelLCBOPIIN::KernelLCBOPIIN() {}

KernelLCBOPIIN::KernelLCBOPIIN(Molecule& _mol): KernelBasic(_mol) {}

KernelLCBOPIIN::~KernelLCBOPIIN() {}

void KernelLCBOPIIN::init() { 
  pForceCalc = new ForceLCBOPIIN(molecule()); 
  pForceCalc->init();
}

void KernelLCBOPIIN::final() { 
  pForceCalc->final();
  delete pForceCalc; 
}

void KernelLCBOPIIN::update() {
  pForceCalc->update();
}

void KernelLCBOPIIN::update(const double* _x) {
  pForceCalc->update(_x);
}

double KernelLCBOPIIN::energy() {
  Molecule& mol = molecule();
  pForceCalc->computeEnergy();
  mol.work(0.0);
  return mol.totalEnergy();
}

void KernelLCBOPIIN::denergy() {//numerical scheme
  pForceCalc->computeForce(1.0e-6);
}

KernelLCBOPIIN::MolTopoPtrs KernelLCBOPIIN::molTopoPtrs() {
  MolTopoPtrs topos; topos.resize(3);
  topos[0] = &pForceCalc->molTopoSR();
  topos[1] = &pForceCalc->molTopoMR();
  topos[2] = &pForceCalc->molTopoLR();
  return topos;
}

/************************************************************/
ForceLCBOPIIN::ForceLCBOPIIN() { reset(); }

ForceLCBOPIIN::ForceLCBOPIIN(AMod::Molecule& _mol): AMod::MolServ(_mol) { 
  reset(); 
}

ForceLCBOPIIN::~ForceLCBOPIIN() {}

void ForceLCBOPIIN::reset() {
  energy.reset();
  topoMR.reset();
  topoLR.reset();
  bSRN1.clear();
  bSRN2.clear();
  aSRN1.clear();
  aMRN2.clear();
  bSRN1IDBak.clear();
  bSRN2IDBak.clear();
  aMRN2IDBak.clear();
  bLRIDBak.clear();
  bSRN1Bak.clear();
  bSRN2Bak.clear();
  aMRN2Bak.clear();
  bLRBak.clear();
}

void ForceLCBOPIIN::init() {
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
  update();
  computeEnergy();
}

void ForceLCBOPIIN::final() {
  reset();
}

void ForceLCBOPIIN::update() {
  AMod::Molecule& mol = molecule();
  mol.topo().update(AMod::CASE_AXES);
  mol.topo().update(AMod::CASE_GENERAL);
  topoMR.update(AMod::CASE_AXES);
  topoMR.update(AMod::CASE_GENERAL);
  topoLR.update(AMod::CASE_AXES);
  topoLR.update(AMod::CASE_GENERAL);
}

void ForceLCBOPIIN::update(const double* _x) {
  AMod::Molecule& mol = molecule();
  int k, kmax = mol.natoms();
  for(k = 0; k < kmax; k++) {
    mol.moveAtom(k,_x+3*k);
    topoMR.moveAtom(k);
    topoLR.moveAtom(k);
  }
}

double ForceLCBOPIIN::computeEnergy() {
  AMod::Molecule& mol = molecule();
  mol.potential(energy.computeE(mol.topo(),topoMR,topoLR));
  return mol.potential();
}

void ForceLCBOPIIN::computeForce(const double& _step) {
  AMod::Molecule& mol = molecule();
  int i, j, natoms = mol.natoms();
  double tmp_x, new_x[3];
  double de, force;
  for(i = 0; i < natoms; i++) {
    atomID = i;
    getNeighbors();
    bak();
    for(j = 0; j < 3; j++) 
      new_x[j] = mol.atom(i)->x[j];
    for(j = 0; j < 3; j++) {
      tmp_x = new_x[j];
      new_x[j] = tmp_x-_step;
      dEnergyMoveAtom(new_x);
      //+step
      new_x[j] = tmp_x+_step;
      de = dEnergyMoveAtom(new_x);
      //force
      force = -de/2.0/_step;
      mol.atom(i)->fpot[j] = force;
      new_x[j] = tmp_x;
    }//end for(j...
    restore();
  }//end for(i...
}

void ForceLCBOPIIN::getNeighbors() {
  int i, imax, j, jmax, k, kmax;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(atomID);
  AMod::Atom& aMR = topoMR.atom(atomID);

  bSRN1.clear(); bSRN2.clear(); 
  aSRN1.clear(); aMRN2.clear();

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

void ForceLCBOPIIN::bak() {
  int k, id;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(atomID);
  AMod::Atom& aMR = topoMR.atom(atomID);
  AMod::Atom& aLR = topoLR.atom(atomID);
  int kmax1 = bSRN1.size();
  int kmax2 = bSRN2.size();
  int kmax3 = aLR.arrows.size();
  int kmax4 = aMRN2.size();

  for(k = 0; k < 3; k++) xBak[k] = aSR->x[k];

  modeBak = mol.topo().mode();
  mol.topo().mode(AMod::MODE_FREEZE);
  topoMR.mode(AMod::MODE_FREEZE);
  topoLR.mode(AMod::MODE_FREEZE);

  bSRN1IDBak.resize(kmax1);
  bSRN1Bak.resize(kmax1);
  for(k = 0; k < kmax1; k++) {
    id = bSRN1[k]->id();
    bSRN1IDBak[k] = id;
    bSRN1Bak[k] = energy.bondDataSR[id];
  }

  bSRN2IDBak.resize(kmax2);
  bSRN2Bak.resize(kmax2);
  for(k = 0; k < kmax2; k++) {
    id = bSRN2[k]->id();
    bSRN2IDBak[k] = id;
    bSRN2Bak[k] = energy.bondDataSR[id];
  }

  bLRIDBak.resize(kmax3);
  bLRBak.resize(kmax3);
  for(k = 0; k < kmax3; k++) {
    id = aLR.arrows[k].bond().id();
    bLRIDBak[k] = id;
    bLRBak[k] = energy.bondDataLR[id];
  }

  aMRN2IDBak.resize(kmax4);
  aMRN2Bak.resize(kmax4);
  for(k = 0; k < kmax4; k++) {
    aMRN2IDBak[k] = aMRN2[k]->id();
    aMRN2Bak[k] = energy.atomDataMR[aMRN2[k]->id()];
  }

  aMRBak = energy.atomDataMR[aMR.id()];
}

void ForceLCBOPIIN::restore() {
  int k, kmax;
  AMod::Molecule& mol = molecule();

  mol.moveAtom(atomID, xBak);
  topoMR.moveAtom(atomID);
  topoLR.moveAtom(atomID);

  mol.topo().mode(modeBak);
  topoMR.mode(modeBak);
  topoLR.mode(modeBak);

  for(k = 0, kmax = bSRN1IDBak.size(); k < kmax; k++)
    energy.bondDataSR[bSRN1IDBak[k]] = bSRN1Bak[k];

  for(k = 0, kmax = bSRN2IDBak.size(); k < kmax; k++)
    energy.bondDataSR[bSRN2IDBak[k]] = bSRN2Bak[k];

  for(k = 0, kmax = bLRIDBak.size(); k < kmax; k++)
    energy.bondDataLR[bLRIDBak[k]] = bLRBak[k];

  for(k = 0, kmax = aMRN2IDBak.size(); k < kmax; k++)
    energy.atomDataMR[aMRN2IDBak[k]] = aMRN2Bak[k];

  energy.atomDataMR[atomID] = aMRBak;
}

double ForceLCBOPIIN::dEnergyMoveAtom(const double* _x) {
  double dE;
  int k, kmax0, kmax1, kmax2, kmax3, kmax4;
  AMod::Molecule& mol = molecule();
  AMod::Atom& aSR = mol.atom(atomID);
  AMod::Atom& aLR = topoLR.atom(atomID);
  AMod::Atom& aMR = topoMR.atom(atomID);

  //old energy for bonds/atoms
  kmax0 = aSRN1.size();
  kmax1 = bSRN1.size();
  kmax2 = bSRN2.size();
  kmax3 = aLR.arrows.size();
  kmax4 = aMRN2.size();

  dE = 0.0;
  for(k = 0; k < kmax1; k++) 
    dE -= energy.bondDataSR[bSRN1[k]->id()].E;
  for(k = 0; k < kmax2; k++) 
    dE -= energy.bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax3; k++) 
    dE -= energy.bondDataLR[aLR.arrows[k].bond().id()].E;
  for(k = 0; k < kmax4; k++)
    dE -= energy.atomDataMR[aMRN2[k]->id()].E;
  dE -= energy.atomDataMR[aMR.id()].E;

  //move atoms
  mol.moveAtom(atomID,_x);
  topoMR.moveAtom(atomID);
  topoLR.moveAtom(atomID);

  //new energy for new/affected bonds/atoms
  for(k = 0; k < kmax1; k++) energy.compute_r0(*bSRN1[k]);
  for(k = 0; k < kmax0; k++) energy.compute_r1(*aSRN1[k]);
  energy.compute_r1(aSR);
  for(k = 0; k < kmax1; k++) energy.compute_r2(*bSRN1[k]);
  for(k = 0; k < kmax2; k++) energy.compute_r2(*bSRN2[k]);
  for(k = 0; k < kmax3; k++) energy.compute_r3(aLR.arrows[k].bond());
  for(k = 0; k < kmax4; k++) energy.compute_r4(*aMRN2[k],mol.topo());
  energy.compute_r4(aMR,mol.topo());

  //sum bonds/atoms energy
  for(k = 0; k < kmax1; k++) 
    dE += energy.bondDataSR[bSRN1[k]->id()].E;
  for(k = 0; k < kmax2; k++) 
    dE += energy.bondDataSR[bSRN2[k]->id()].E;
  for(k = 0; k < kmax3; k++) 
    dE += energy.bondDataLR[aLR.arrows[k].bond().id()].E;
  for(k = 0; k < kmax4; k++)
    dE += energy.atomDataMR[aMRN2[k]->id()].E;
  dE += energy.atomDataMR[aMR.id()].E;

  return dE;
}

}/* AMod */
