#include "Energy.h"
#include "S.h"
#include "handy.h"
#include "globalConstants.h"
#include <cmath>

namespace LCBOPII {

static double fVmr0(double r);
static double fVmr1(double r);
static double fVmr2(double r);

double Energy::fSVmr(const AMod::Arrow& arrowMR,
		     const AMod::MolTopo& topoSR) {
  int i = arrowMR.host().id();
  int j = arrowMR.moon().id();
  const AMod::Atom& aiSR = topoSR.atom(i);
  const AMod::Arrow& arijMR = arrowMR;
  int aiSR_n = aiSR.arrows.size();
  /* Nij & xdbij */
  int k, m, akSR_n;
  double Suelki, Nelki;
  double xdbij, Ndbij, Nki, Mki, Nij;
  Nij = 0.0;
  Ndbij = 4.0;
  for(k = 0; k < aiSR_n; k++) {
    const AMod::Arrow& arikSR = aiSR.arrows[k];
    const AMod::Atom& akSR = arikSR.moon();
    if(akSR.id() == j) continue;
    /* Mki */
    Mki = 0.0;
    akSR_n = akSR.arrows.size();
    for(m = 0; m < akSR_n; m++) {
      const AMod::Arrow& arkmSR = akSR.arrows[m];
      if(arkmSR.moon().id() == i) continue;
      const BondDatumSR& bdkm = bondDataSR[arkmSR.bond().id()];
      Mki += bdkm.SdN*bdkm.M_[1-arkmSR.position()];
    }//end for(m...
    /* Nelki */
    const BondDatumSR& bdik = bondDataSR[arikSR.bond().id()];
    Nki = bdik.N_[1-arikSR.position()];
    Suelki = bdik.Suel_[1-arikSR.position()];
    Nelki = (1.0-Suelki)*(4.0-Mki)/(Nki+1.0-Mki)+Suelki*4.0/(Nki+1.0);
    if(Nelki >= 3.5) Nelki = 3.0;
    else if(Nelki >= 2.5) Nelki = 3.0-0.5*(3.5-Nelki)*(3.5-Nelki);
    /* Ndbij */
    Ndbij -= bdik.SdN*Nelki;
    /* Nij */
    Nij += bdik.SdN;
  }//end for(k...
  if(Ndbij>3) return 0.0;
  if(Ndbij<0) xdbij = 0.0;
  else xdbij = (Ndbij-int(Ndbij));
  /* gammaij */
  double NijP, gammaij, sum, temp, temp2;
  for(k = 0, sum = 0.0; k < aiSR_n; k++) {
    const AMod::Arrow& arikSR = aiSR.arrows[k];
    if(arikSR.moon().id() == j) continue;
    temp = 1.0+costh(arijMR.dx,arikSR.dx);
    temp2 = temp*temp;
    sum += temp2*temp2*bondDataSR[arikSR.bond().id()].SdN;
  }
  NijP = (Nij<1 ? 0.5*(1.0+Nij*Nij) : Nij);
  gammaij = 1.0/(1.0+(Const::B/NijP)*sum);
  double gamma0 = fSugamma0(gammaij);
  double gamma2 = fSugamma2(gammaij);
  if(gamma2 == 0.0) return 0.0;
  /* Vmrij */
  double Vmrij, rij;
  rij = arijMR.bond().r;
  temp = fSddb(xdbij);
  if(Ndbij < 1.0) {
    Vmrij = 
      temp*gamma0*fVmr0(rij)+(1.0-temp)*square(gamma2)*fVmr1(rij);
  }
  else if(Ndbij < 2.0) {
    Vmrij = 
      temp*square(gamma2)*fVmr1(rij)+
      (1.0-temp)*gamma2*fVmr2(rij);
  }
  else {
    Vmrij = temp*gamma2*fVmr2(rij);
  }
  return fSumr(rij)*Vmrij;
}

double fVmr0(double r) {
  double dr = Const::rmr1-r;
  if(dr <= 0.0) return 0.0;
  return Const::Amr0*dr*dr*dr;
}

double fVmr1(double r) {
  double dr = Const::rmr1-r;
  if(dr <= 0.0) return 0.0;
  return Const::Amr1*dr*dr*dr;
}

double fVmr2(double r) {
  double dr = Const::rmr2-r;
  if(dr <= 0.0) return 0.0;
  return Const::Amr2*dr*dr;
}

}/* LCBOPII */
