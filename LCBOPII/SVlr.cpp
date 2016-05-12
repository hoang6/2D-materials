#include "Energy.h"
#include "S.h"
#include "globalConstants.h"
#include <cmath>

namespace LCBOPII {

static double fVlr1(double r);
static double fVlr2(double r);

double Energy::fSVlr(const AMod::Bond& bondLR) {
  if(bondLR.r < Const::r0) 
    return (1.0-fSdsr(bondLR.r))*fVlr1(bondLR.r)*fSdlr(bondLR.r);
  return (1.0-fSdsr(bondLR.r))*fVlr2(bondLR.r)*fSdlr(bondLR.r);
}

double fVlr1(double r) {
  return 
    Const::eps1*(std::exp(-2.0*Const::lambda1*(r-Const::r0))-
		 2.0*std::exp(-Const::lambda1*(r-Const::r0)))+Const::v1;
}

double fVlr2(double r) {
  return 
    Const::eps2*(std::exp(-2.0*Const::lambda2*(r-Const::r0))-
		 2.0*std::exp(-Const::lambda2*(r-Const::r0)))+Const::v2;
}

}/* LCBOPII */
