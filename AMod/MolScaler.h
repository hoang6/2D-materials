#ifndef _AMOD_MOL_SCALER_H
#define _AMOD_MOL_SCALER_H

#include "MolService.h"

namespace AMod {

class MolScaler: public ConstMolServ {
public:
  MolScaler();
  MolScaler(Molecule& _mol);
  void reset();
  virtual void attach(const Molecule& _mol);
  void init();
  void scale(int _atomID, double _sx[3]) const;

protected:
  double h[9], inv_h[9];
};

}/* AMod */

#endif
