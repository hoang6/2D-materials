#ifndef _MD_DATA_H
#define _MD_DATA_H

#include "../AMod/Molecule.h"
#include "../AMod/KernelBasic.h"
#include "../AMod/PTE.h"
#include <vector>

namespace MD {

class Data {
public:
  AMod::PotType pot;
  AMod::Molecule mol;
  typedef std::vector<double> Array; Array vel;
  //------------------------------------------------------------
  Data() { reset(); }
  Data(AMod::PotType _pot, AMod::Molecule& _mol, Array& _vel): 
    pot(_pot), mol(_mol), vel(_vel) {}
  void reset() { pot = AMod::REBO; mol.reset(); vel.clear(); }
  double mass(int _atomID) const {
    return AMod::PTE::atomicMass(mol.atom(_atomID)->type)*103.642722245; }
  void forceConsistency() { vel.resize(3*mol.natoms(),0.0); }
  void resetTotalMomentum();
  void rescaleVel(double _kT);
  void rmFreeC(const std::vector<int>& _freeAtomIDs);

private:
  struct Point {double x[3]; Util::PtrCont<Point>::ID id; };
};

}/* MD */


#endif
