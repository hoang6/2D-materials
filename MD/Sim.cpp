#include "Sim.h"

namespace MD {

Sim::Sim() {}

Sim::~Sim() {}

int Sim::run(NVE& _ensem, IO& _io) {
  int i, step = 0;
  if(!_io.findKey("DTIME")) return 1;
  if(!_io.findKey("MAX_STEP")) return 2;
  if(!_io.findKey("DUMP_STEP")) return 3;
  if(_io.openFiles() != 0) return 4;
  for(step = 1; step <= _io.maxStep; step++) {
    _ensem.step(_io.dtime);
    if(step%_io.dumpStep == 0) {
      _ensem.molecule().io().dumpTxt(_io.foutMol); 
      _io.foutMol << std::endl;
      for(i = 0; i < _ensem.ndata(); i++)
	_io.foutVel << _ensem.vel(i) << " ";
      _io.foutVel << std::endl;
      _ensem.dumpEnergy(_io.foutEne);
    }
    if(this->stop(_ensem,step)) return 0;
  }
  _io.closeFiles();
  getFreeAtomIDs(_ensem);
  return 0;
}

bool Sim::stop(NVE& _ensem, int _step) { return false; }

void Sim::getFreeAtomIDs(NVE& _ensem) {
  Data& data = _ensem.data();
  AMod::KernelBasic::MolTopoPtrs topoPtrs = _ensem.kernel().molTopoPtrs();
  fAtomIDs.clear();
  for(int i = 0, imax = data.mol.natoms(); i < imax; i++) {
    int nNeighbs = 0;
    for(int j = 0, jmax = topoPtrs.size(); j < jmax; j++)
      nNeighbs += topoPtrs[j]->atom(i).arrows.size();
    if(nNeighbs == 0) fAtomIDs.push_back(i);
  }
}

}/* MD */
