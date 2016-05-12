#ifndef _MC_MC_AGENT_CTRLZ_H
#define _MC_MC_AGENT_CTRLZ_H

#include "../MC/MC.h"

namespace MC {

/** "CtrlZ" stands for "control z" **/
class MCAgent_CtrlZ: public MCAgent {
public:
  MCAgent_CtrlZ();
  MCAgent_CtrlZ(AMod::Molecule& _mol);
  virtual ~MCAgent_CtrlZ();
  virtual void pertAtom(double _dxmax, int _id, double _new_x[3]);

  Util::Accessors<double> zmin;
  Util::Accessors<double> zmax;
};

}/* MC */

#endif
