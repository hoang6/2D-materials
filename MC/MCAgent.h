#ifndef _MC_MC_AGENT_H
#define _MC_MC_AGENT_H

#include "../AMod/MolService.h"
#include "../Util/Accessors.h"
#include "../Util/RNG.h"
#include "MCBasic.h"

namespace MC {

class MCAgent: public AMod::MolServ {
public:
  MCAgent();
  MCAgent(AMod::Molecule& _mol);
  virtual ~MCAgent();
  virtual void init();
  virtual void final();
  virtual void pertAtom(double _dxmax, int _id, double _new_x[3]);
  virtual void accPertAtom(bool _flag);
  virtual void pertAxes(double _dhmax);
  virtual void accPertAxes(bool _flag);

protected:
  double urand(double _a, double _b) const;
};

}/* MC */

#endif
