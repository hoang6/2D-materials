#ifndef _SW_SW_BARRIER_H
#define _SW_SW_BARRIER_H

#include "SWStimulus.h"
#include "../AMod/MolService.h"

namespace SW {

class SWBarrier: public AMod::MolServ {
public:
  SWBarrier();
  void init(KernelType _kernelType, const SWCreator::KeyAtoms& _keyAtoms);
  void final();
  void optimize(double _reactionCoord,
		AMod::MolTopoMode _mode = AMod::MODE_FREEZE);
};

}/* SW */

#endif
