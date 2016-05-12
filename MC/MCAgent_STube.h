#ifndef _MC_MC_AGENT_STUBE_H
#define _MC_MC_AGENT_STUBE_H

#include "../MC/MC.h"
#include <vector>

namespace MC {

/** "STube" stands for "smooth tube"
    Note no atom should be fixed **/
class MCAgent_STube: public MCAgent {
public:
  MCAgent_STube();
  MCAgent_STube(AMod::Molecule& _mol);
  virtual ~MCAgent_STube();
  virtual void init();
  virtual void pertAtom(double _dxmax, int _id, double _new_x[3]);
  virtual void accPertAtom(bool _flag);
  virtual void pertAxes(double _dhmax);
  virtual void accPertAxes(bool _flag);

private:
  double _r, _lz;
  std::vector<double> _th;
  int _thFixID;

  double _new_r, _new_lz;
  double _new_id, _new_th;
};

}/* MC */

#endif
