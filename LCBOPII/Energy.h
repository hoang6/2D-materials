#ifndef _LCBOPII_ENERGY_H
#define _LCBOPII_ENERGY_H

#include "../AMod/Molecule.h"
#include "../Util/PtrCont.h"
#include <vector>

namespace LCBOPII {

class Energy;
class BondDatumSR;
class BondDatumLR; 
class AtomDatumMR;
typedef Util::PtrCont<BondDatumSR> BondDataSR;
typedef Util::PtrCont<BondDatumLR> BondDataLR;
typedef Util::PtrCont<AtomDatumMR> AtomDataMR;

class BondDatumSR {
public:
  double SdN;
  bool fullNeighb; //SdN == 1
  double N_[2];    //N_[*] = N[*]-SdN
  double M_[2];    //M_[*] = fSuM(N_[*])
  double Suel_[2];
  double E;
  Util::PtrCont<BondDatumSR>::ID id;
};

class BondDatumLR {
public:
  double E;
  Util::PtrCont<BondDatumLR>::ID id;
};

class AtomDatumMR {
public:
  double E;
  Util::PtrCont<AtomDatumMR>::ID id;
};

class Energy {
public:
  BondDataSR bondDataSR;
  BondDataLR bondDataLR;
  AtomDataMR atomDataMR;

  Energy();
  Energy(const Energy& energy);
  ~Energy();
  Energy& operator= (const Energy& energy);
  void reset();
  static AMod::BondRange bondRangeSR();
  static AMod::BondRange bondRangeMR();
  static AMod::BondRange bondRangeLR();
  static double cellLengthSR();
  static double cellLengthMR();
  static double cellLengthLR();
  void compute_r0(const AMod::Bond& bondSR);    //SdN
  void compute_r1(const AMod::Atom& atomSR);    //N_
  void compute_r2(const AMod::Bond& bondSR);    //SVsr
  void compute_r3(const AMod::Bond& bondLR);    //SVlr
  void compute_r4(const AMod::Atom& atomMR,
		  const AMod::MolTopo& topoSR); //Emr
  double computeE(const AMod::MolTopo& topoSR,
		  const AMod::MolTopo& topoMR,
		  const AMod::MolTopo& topoLR);
  /********* valid only after computeE(...) *********/
  double getESR() const { return potSR; }
  double getEMR() const { return potMR; }
  double getELR() const { return potLR; }
  double getEtot() const { return potSR+potMR+potLR; }

private:
  struct FDatum {double W; int tN; double tM; std::vector<int> ccArIDs;};
  void fPreF(const AMod::Arrow& arrowSR, std::vector<FDatum>& fdata);
  double fFAT(const AMod::Bond& bondSR);
  double fb(const AMod::Arrow& arrowSR);
  double fSVsr(const AMod::Bond& bondSR);
  double fSVlr(const AMod::Bond& bondLR);
  double fSVmr(const AMod::Arrow& arrowMR, 
	       const AMod::MolTopo& topoSR);

  double potSR, potMR, potLR;
};

}/* LCBOP II */

#endif
