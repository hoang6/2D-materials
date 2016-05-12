#include "MCAgent_CtrlZ.h"

namespace MC {

MCAgent_CtrlZ::MCAgent_CtrlZ() {}

MCAgent_CtrlZ::MCAgent_CtrlZ(AMod::Molecule& _mol): MCAgent(_mol) {}

MCAgent_CtrlZ::~MCAgent_CtrlZ() {}

void MCAgent_CtrlZ::pertAtom(double _dxmax, int _id, double _new_x[3]) {
  double atomz, zbmin, zbmax;
  const AMod::Atom& atom = molecule().atom(_id);
  for(int k = 0; k < 3; k++) {
    if(atom->fixed[k]) 
      _new_x[k] = atom->x[k];
    else if(k != 2)
      _new_x[k] = atom->x[k]+_dxmax*urand(-1,1);
    else {
      atomz = atom->x[k];
      zbmin = atomz-_dxmax;
      zbmax = atomz+_dxmax;
      if(zbmin < zmin()) zbmin = zmin();
      if(zbmax > zmax()) zbmax = zmax();
      _new_x[k] = urand(zbmin,zbmax);
    }
  }//end for(int k...
}

}/* MC */
