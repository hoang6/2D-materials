#ifndef _AMOD_MOL_FACTORY_H
#define _AMOD_MOL_FACTORY_H

#include "MolService.h"
#include "Box.h"
#include <vector>
#include <utility>

namespace AMod {

class MolFactory: public MolServ {
public:
  MolFactory();
  MolFactory(Molecule& _mol);
  void reset();
  /****** factory ******/
  void makeGPH_PBC(const int chiV[2], const int heightV[2]);
  void makeGPH_NPBC(const int chiV[2], const int heightV[2], bool _capH = true);
  void makeCNT_PBC(const int chiV[2], const int heightV[2]);
  void makeCNT_NPBC(const int chiV[2], const int heightV[2], bool _capH = true);
  void makeC60();
  /****** box edges ******/
  const Box::Information& boxInformation() const;

protected:
  void makeGPH(const int chiV[2], const int heightV[2]);
  void makeCNT(const int chiV[2], const int heightV[2]);
  void clearBoundaryBonds() const;
  void capHydrogen() const;
  double findHeight(int chic, int chir, double height);
  static void rollup(const double& r, double* xyz);

  Box::Information boxInfo;
};

/************************************************************/
inline const Box::Information& MolFactory::boxInformation() const { return boxInfo; }

} /* AMod */

#endif
