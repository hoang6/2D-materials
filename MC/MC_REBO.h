#ifndef _MC_MC_REBO_H
#define _MC_MC_REBO_H

#include "MCBasic.h"
#include "../REBO/Energy.h"
#include <vector>
#include <set>

namespace MC {

class MC_REBO: public MCBasic { 
public:
  MC_REBO();
  MC_REBO(AMod::Molecule& _mol);
  virtual ~MC_REBO();
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
  AMod::MolTopo& molTopo();
  const AMod::MolTopo& molTopo() const;

  Util::Accessors<REBO::Energy> energy;

private:
  void clear();
  void bakAffected();
  void bakDeleted();
  void restoreAffected();
  void restoreDeleted();
  void getNeighbors(int _atomID);

  /********** for getNeighbBonds(int) **********/
  typedef std::vector<AMod::Bond*> BondPtrs;
  typedef std::vector<AMod::Atom*> AtomPtrs;
  BondPtrs bN1;
  BondPtrs bN2;
  /********** for bak **********/
  AMod::Molecule molBak;
  REBO::Energy energyBak;
  double potentialBak;
  int atomIDBak;
  double xBak[3];
  AMod::Atom::Data atomDataBak;
  std::vector<int> bondDataIDBak;
  std::vector<REBO::BondDatum> bondDataBak;
  std::vector<int> bondBak;
  std::vector<REBO::BondDatum> bondDataDelBak;
};

}/* MC */

#endif
