#include "Energy.h"
#include "handy.h"
#include "globalConstants.h"
#include <cmath>

namespace LCBOPII {

namespace {

/********** H **********/
double fH(double x);
/********** G **********/
double fy0(double z);
double fG1(double y);
double dfG1(double y);
double fG2(double y, double z, double y0);
double fG(double y, double z);

}/* unnamed */

double Energy::fb(const AMod::Arrow& arrowSR) {
  const AMod::Atom& ai = arrowSR.host();
  const AMod::Arrow& arij = arrowSR;
  const BondDatumSR& bdij = bondDataSR[arij.bond().id()];
  const double rij = arij.bond().r;
  const double Nij = bdij.N_[arij.position()];
  double costhijk;
  double Nijk;
  int k, kmax = ai.arrows.size();
  double sum = 0.0;
  for(k = 0; k < kmax; k++) {
    if(k == arij.id()) continue;
    const AMod::Arrow& arik = ai.arrows[k];
    const BondDatumSR& bdik = bondDataSR[arik.bond().id()];
    costhijk = costh(arij.dx,arik.dx);
    Nijk = Nij-bdik.SdN;
    sum += bdik.SdN*fH(rij-arik.bond().r)*fG(costhijk,Nijk);
  }
  return 1.0/std::sqrt(1.0+sum);
}

namespace {

/********** H **********/
double fH(double x) {
  if(x < -Const::d) {
    double temp = Const::kappa*(x+Const::d);
    temp = square(square(temp));
    return Const::L*(1+Const::kappa*(x+Const::d)*std::pow(1.0+temp,-0.25));
  }
  if(x <= Const::d) {
    double x2 = x*x;
    double x4 = x2*x2;
    double x6 = x2*x4;
    return 1.0+Const::C1*x+0.5*Const::C12*x2+Const::C4*x4+Const::C6*x6;
  }
  return Const::R0+Const::R1*(x-Const::d);
}

/*************** G ***************/
double fy0(double z) {
  return Const::Ay0+Const::By0*(z+z*z);
}

double fG1(double y) {
  double sum = 0.0;
  double ypow = 1.0;
  if(y < -0.5) {
    for(int k = 0; k < 3; k++) {
      sum += Const::g1[k]*ypow;
      ypow *= y;
    }
    return Const::gmin+(y+1.0)*(y+1.0)*sum;
  }
  if(y < -1.0/3.0) {
    for(int k = 0; k < 5; k++) {
      sum += Const::g2[k]*ypow;
      ypow *= y;
    }
    return Const::ggr+(y+0.5)*sum;
  }
  /* -1.0/3.0 <= y && y <= 1.0 */
  for(int k = 0; k < 5; k++) {
    sum += Const::g3[k]*ypow;
    ypow *= y;
  }
  return Const::gmax+(y-1.0)*(y-1.0)*sum;
}

double dfG1(double y) {
  double sum1 = 0.0; 
  double sum2 = 0.0;
  double ypow = 1.0;
  if(y < -0.5) {
    for(int k = 0; k < 3; k++) {
      sum1 += Const::g1[k]*ypow;
      if(k < 2) sum2 += Const::g1[k+1]*(k+1)*ypow;
      ypow *= y;
    }
    return 2.0*(y+1.0)*sum1+(y+1.0)*(y+1.0)*sum2;
  }
  if(y < -1.0/3.0) {
    for(int k = 0; k < 5; k++) {
      sum1 += Const::g2[k]*ypow;
      if(k < 4) sum2 += Const::g2[k+1]*(k+1)*ypow;
      ypow *= y;
    }
    return sum1+(y+0.5)*sum2;
  }
  /* -1.0/3.0 <= y && y <= 1.0 */
  for(int k = 0; k < 5; k++) {
    sum1 += Const::g3[k]*ypow;
    if(k < 4) sum2 += Const::g3[k+1]*(k+1)*ypow;
    ypow *= y;
  }
  return 2.0*(y-1.0)*sum1+(y-1.0)*(y-1.0)*sum2;
}

double fG2(double y, double z, double y0) {
  double z4 = z*z*z*z;
  double one_y02 = (1.0-y0)*(1.0-y0);
  double one_y03 = one_y02*(1.0-y0);
  double G1y0 = fG1(y0);

  double gzmax = Const::gmax-(Const::Ag+Const::Bg*z+Const::Cg*z4)*one_y03;
  double gz[3];
  gz[2] = Const::Dg*z4/(1.0+Const::Eg*z4);
  gz[1] = dfG1(y0)/one_y02+2.0*(G1y0-gzmax)/one_y03-2.0*gz[2]*y0;
  gz[0] = (G1y0-gzmax)/one_y02-gz[1]*y0-gz[2]*y0*y0;

  return gzmax+(1.0-y)*(1.0-y)*(gz[0]+gz[1]*y+gz[2]*y*y);
}

double fG(double y, double z) {
  double y0 = fy0(z);
  if(y0 > y) return fG1(y);
  return fG2(y,z,y0);
}

}/* unnamed */

}/* LCBOPII */
