#ifndef _AMOD_KERNEL_LCBOPIIN_H
#define _AMOD_KERNEL_LCBOPIIN_H

#include "KernelBasic.h"
#include "../LCBOPII/Energy.h"

namespace AMod {

class ForceLCBOPIIN;

class KernelLCBOPIIN: public KernelBasic {
public:
  KernelLCBOPIIN();
  KernelLCBOPIIN(Molecule& _mol);
  virtual ~KernelLCBOPIIN();
  virtual void init();
  virtual void final();
  virtual void update();
  virtual void update(const double* _x);
  virtual double energy();
  virtual void denergy();
  virtual MolTopoPtrs molTopoPtrs();
  ForceLCBOPIIN& kernelF() { return *pForceCalc; }

protected:
  ForceLCBOPIIN* pForceCalc;
};

/************************************************************/
class ForceLCBOPIIN: public AMod::MolServ {
public:
  ForceLCBOPIIN();
  ForceLCBOPIIN(AMod::Molecule& _mol);
  ~ForceLCBOPIIN();
  void reset();
  void init();
  void final();
  void update();
  void update(const double* _x);
  double computeEnergy();
  void computeForce(const double& _step); //only called after computeEnergy(...)
  //------------------------------------------------------------
  LCBOPII::Energy& kernelE() { return energy; }
  AMod::MolTopo& molTopoSR() { return molecule().topo(); }
  AMod::MolTopo& molTopoMR() { return topoMR; }
  AMod::MolTopo& molTopoLR() { return topoLR; }

protected:
  LCBOPII::Energy energy;
  //topoSR is molecule().topo()
  AMod::MolTopo topoMR;
  AMod::MolTopo topoLR;

private:
  void getNeighbors();
  void bak();
  void restore();
  double dEnergyMoveAtom(const double* _dx);

  int atomID;
  /********** for getNeighbors(int) **********/
  typedef std::vector<AMod::Bond*> BondPtrs;
  typedef std::vector<AMod::Atom*> AtomPtrs;
  BondPtrs bSRN1; //1st neighbor SR bonds
  BondPtrs bSRN2; //affected SR bonds
  AtomPtrs aSRN1; //1st neighbor SR atoms
  AtomPtrs aMRN2; //affected MR atoms
  /********** for bak **********/
  double xBak[3];
  AMod::MolTopoMode modeBak;
  std::vector<int> bSRN1IDBak;
  std::vector<int> bSRN2IDBak;
  std::vector<int> aMRN2IDBak;
  std::vector<int> bLRIDBak;
  std::vector<LCBOPII::BondDatumSR> bSRN1Bak;
  std::vector<LCBOPII::BondDatumSR> bSRN2Bak;
  std::vector<LCBOPII::AtomDatumMR> aMRN2Bak;
  std::vector<LCBOPII::BondDatumLR> bLRBak;
  LCBOPII::AtomDatumMR aMRBak;
};

}/* AMod */

#endif
