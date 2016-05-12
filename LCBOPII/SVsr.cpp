#include "Energy.h"
#include "S.h"
#include "globalConstants.h"
#include <cmath>

namespace LCBOPII {

double Energy::fSVsr(const AMod::Bond& bondSR) {
  double VsrR, VsrA;
  double b[2], FAT, B;
  VsrR = Const::Asr*std::exp(-Const::alpha*bondSR.r);
  VsrA = 
    Const::Bsr1*std::exp(-Const::beta1*bondSR.r)+
    Const::Bsr2*std::exp(-Const::beta2*bondSR.r);
  b[0] = fb(bondSR.arrow(0));
  b[1] = fb(bondSR.arrow(1));
  FAT = fFAT(bondSR);
  B = 0.5*(b[0]+b[1])+FAT;
  return fSdsr(bondSR.r)*(VsrR-B*VsrA);
}

}/* LCBOPII */
