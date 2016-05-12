#ifndef _MC_MC_BASIC_H
#define _MC_MC_BASIC_H

#include "../AMod/AMod.h"

namespace MC {

class MCBasic: public AMod::MolServ {
public:
  MCBasic() {};
  MCBasic(AMod::Molecule& _mol): AMod::MolServ(_mol) {}
  virtual ~MCBasic() {};
  virtual void init() = 0;
  virtual void final() = 0;
  virtual void save() = 0;
  virtual void restore() = 0;  
  virtual double computeEnergy() = 0;
  virtual double dEnergyNewAtom(int _atomID, 
				const AMod::Atom::Data& _atomData) = 0;
  virtual double dEnergyDeleteAtom(int _atomID) = 0;
  virtual double dEnergyMoveAtom(int _atomID, const double* _x) = 0;
  virtual void unNewAtom() = 0;
  virtual void unDeleteAtom() = 0;
  virtual void unMoveAtom() = 0;
};

}/* MC */

#endif
