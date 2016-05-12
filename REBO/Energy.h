#ifndef _REBO_ENERGY_H
#define _REBO_ENERGY_H

#include "../AMod/Molecule.h"
#include "../Util/PtrCont.h"
#include <iostream>

namespace REBO {

class Energy;
class BondDatum;
typedef Util::PtrCont<BondDatum> BondData;

class BondDatum {
public:
  //-------------------- double sided --------------------
  //depend on types
  int type[2];
  double bsigmapi[2];
  double N[2][3], F[2];
  //-------------------- single sided --------------------
  //depend on types
  double Q, A, alpha, B[3], beta[3];
  double sum_theta2;
  double bDH, bbar;
  double fc;
  double sumFfc[2];
  double Nconj, Fij, Tij;
  double VApre;
  double VR, VA, E;
  double f[3];
  Util::PtrCont<BondDatum>::ID id;
}; 

class Energy {
public:
  BondData bondData;

  Energy();
  Energy(const Energy& energy);
  ~Energy();
  Energy& operator= (const Energy& energy);
  void reset();
  static AMod::BondRange bondRange();
  static double cellLength();
  void compute_r0(const AMod::Bond& bond);
  void compute_r1(const AMod::Bond& bond);
  void compute_r2(const AMod::Bond& bond);
  double computeE(const AMod::MolTopo& topo);
  void dumpBondData(std::ostream& fdebug, const AMod::Bond& bond) const;
};

} /* REBO */

#endif
