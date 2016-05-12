#ifndef _MD_SIM_H
#define _MD_SIM_H

#include "IO.h"

namespace MD {

class Sim {
public:
  Sim();
  virtual ~Sim();
  int run(NVE& _ensem, IO& _io);
  virtual bool stop(NVE& _ensem, int _step);
  const std::vector<int>& freeAtomIDs() { return fAtomIDs; }

private:
  void getFreeAtomIDs(NVE& _ensem);
  std::vector<int> fAtomIDs;
};

}/* MD */

#endif
