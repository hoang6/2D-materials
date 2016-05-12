#ifndef _MC_MC_LCBOPII_H
#define _MC_MC_LCBOPII_H

#include "MCBasic.h"
#include "../LCBOPII/Energy.h"
#include <vector>

namespace MC {

class MC_LCBOPII: public MCBasic { 
public:
  MC_LCBOPII();
  MC_LCBOPII(AMod::Molecule& _mol);
  virtual ~MC_LCBOPII();
  virtual void init();
  virtual void final();
  virtual void save();
  virtual void restore();  
  virtual double computeEnergy();
  virtual double dEnergyNewAtom(int _atomID, 
				const AMod::Atom::Data& _atomData);
  virtual double dEnergyDeleteAtom(int _atomID);
  virtual double dEnergyMoveAtom(int _atomID, const double* _x);
  virtual void unNewAtom();
  virtual void unDeleteAtom();
  virtual void unMoveAtom();
  AMod::MolTopo& molTopoSR();
  AMod::MolTopo& molTopoMR();
  AMod::MolTopo& molTopoLR();
  const AMod::MolTopo& molTopoSR() const;
  const AMod::MolTopo& molTopoMR() const;
  const AMod::MolTopo& molTopoLR() const;

  Util::Accessors<LCBOPII::Energy> energy;

protected:
  /********** for LCBOPII **********/
  //topoSR is molecule().topo()
  AMod::MolTopo topoMR;
  AMod::MolTopo topoLR;  

private:
  void clear();
  void bakAffected();
  void bakDeleted();
  void restoreAffected();
  void restoreDeleted();
  void getNeighbors(int _atomID);

  /********** for getNeighbors(int) **********/
  typedef std::vector<AMod::Bond*> BondPtrs;
  typedef std::vector<AMod::Atom*> AtomPtrs;
  BondPtrs bSRN1; //1st neighbor SR bonds
  BondPtrs bSRN2; //affected SR bonds
  AtomPtrs aSRN1; //1st neighbor SR atoms
  AtomPtrs aMRN2; //affected MR atoms
  /********** for bak **********/
  AMod::Molecule molBak; 
  AMod::MolTopo topoMRBak;
  AMod::MolTopo topoLRBak;
  LCBOPII::Energy energyBak;
  double potentialBak;
  int atomIDBak;
  double xBak[3];
  AMod::Atom::Data atomDataBak;
  std::vector<int> bondDataSRIDBak;
  std::vector<int> atomDataMRIDBak;
  std::vector<LCBOPII::BondDatumSR> bondDataSRBak;
  std::vector<LCBOPII::AtomDatumMR> atomDataMRBak;
  std::vector<int> bondSRBak;
  std::vector<int> bondLRBak;
  std::vector<int> bondMRBak;
  std::vector<LCBOPII::BondDatumSR> bondDataSRDelBak;
  std::vector<LCBOPII::BondDatumLR> bondDataLRDelBak;
  LCBOPII::AtomDatumMR atomDatumMRDelBak;
};

}/* MC */

#endif
