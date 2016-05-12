#ifndef _SW_SW_CREATOR_H
#define _SW_SW_CREATOR_H

#include "SWStimulus.h"
#include "../AMod/MolService.h"

namespace SW {

class SWCreator: public AMod::MolServ {
public:
  typedef struct { int id[2][3]; } KeyAtoms;

  SWCreator();
  SWCreator(AMod::Molecule& _mol);
  /****** change topology ******/
  SWStimulus act(const SWStimulus& _stim);
  static KeyAtoms quad(const AMod::Molecule& _mol, const SWStimulus& _stim);
};

} /* SW */

#endif
