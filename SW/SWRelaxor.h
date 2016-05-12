#ifndef _SW_SW_RELAXOR_H
#define _SW_SW_RELAXOR_H

#include "../AMod/Atom.h"
#include "../AMod/MolService.h"
#include "../Util/util.h"
#include <vector>

namespace SW {

class SWRelaxor: public AMod::MolServ {
public:
  SWRelaxor();
  SWRelaxor(AMod::Molecule& _mol);
  double relax(int _bondID, double _tol = 1.0e-9);

protected:
  void setup(int _bondID);
  int nDOF() const;
  double func(double* _x);
  void dfunc(double* _x, double* _g);
  static double stdBondLength(int type1, int type2);
  static double stdBondLength(const AMod::Bond& bond);

  std::vector<int> atomIDs;
  std::vector<AMod::DInt> atomIDPairs;
};

}/* SW */

#endif
