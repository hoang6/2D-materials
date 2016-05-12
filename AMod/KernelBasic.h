#ifndef _AMOD_KERNEL_BASIC_H
#define _AMOD_KERNEL_BASIC_H

#include "MolService.h"
#include "Molecule.h"
#include <iostream>
#include <vector>

namespace AMod {

typedef enum {REBO, LCBOPIIN} PotType;

class KernelBasic: public MolServ {
public:
  KernelBasic();
  KernelBasic(Molecule& _mol);
  virtual ~KernelBasic();
  virtual void init();
  virtual void final();
  virtual void update() = 0; //update topology
  virtual void update(const double* _x) = 0;
  virtual double energy() = 0; //the topology must be correct
  virtual void denergy() = 0; //may depend on energy()
  typedef std::vector<MolTopo*> MolTopoPtrs;
  virtual MolTopoPtrs molTopoPtrs();
};

/************************************************************/
inline KernelBasic::KernelBasic() {}

inline KernelBasic::KernelBasic(Molecule& _mol): MolServ(_mol) {}

inline KernelBasic::~KernelBasic() {}

inline void KernelBasic::init() {}

inline void KernelBasic::final() {}

inline KernelBasic::MolTopoPtrs KernelBasic::molTopoPtrs() {
  return MolTopoPtrs(1,&molecule().topo());
}

}/* AMod */

#endif
