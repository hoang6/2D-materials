#ifndef _AMOD_KERNEL_REBO_H
#define _AMOD_KERNEL_REBO_H

#include "KernelBasic.h"
#include "../REBO/Energy.h"
#include "../REBO/Force.h"

namespace AMod {

class KernelREBO: public KernelBasic {
public:
  KernelREBO();
  KernelREBO(Molecule& _mol);
  virtual ~KernelREBO();
  virtual void init();
  virtual void final();
  virtual void update();
  virtual void update(const double* _x);
  virtual double energy();
  virtual void denergy();

protected:
  REBO::Energy* pEnergyCalc;
  REBO::Force* pForceCalc;
};

}/* AMod */

#endif
