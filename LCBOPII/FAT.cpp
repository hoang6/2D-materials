#include "Energy.h"
#include "S.h"
#include "handy.h"
#include "globalConstants.h"
#include <vector>
#include <cmath>

namespace LCBOPII {

namespace {

double ft(double y, double z, double Delta);
void ftvec(const double rij[3], const double rik1[3], const double rik2[3], double t[3]);

}/* unnamed */

void Energy::fPreF(const AMod::Arrow& arrowSR, std::vector<FDatum>& fdata) {
  const AMod::Arrow& arij = arrowSR;
  const AMod::Atom& ai = arij.host();
  int c, ai_n = ai.arrows.size();

  /* Initialize fullArID & fracArID */
  std::vector<int> fullArIDi;
  std::vector<int> fracArIDi;
  fullArIDi.reserve(ai_n);
  fracArIDi.reserve(ai_n);
  for(c = 0; c < ai_n; c++) {
    if(c == arij.id()) continue;
    if(bondDataSR[ai.arrows[c].bond().id()].fullNeighb)
      fullArIDi.push_back(c);
    else
      fracArIDi.push_back(c);
  }
  int nfullArIDi = fullArIDi.size();
  int nfracArIDi = fracArIDi.size();
  
  /* compute fdata */
  fdata.clear();
  if(nfullArIDi < 3) {
    //init
    fdata.resize(1);
    fdata[0].W = 1.0;
    fdata[0].tN = nfullArIDi;
    fdata[0].tM = 0.0;
    for(c = 0; c < nfullArIDi; c++) {
      const AMod::Arrow& arik = ai.arrows[fullArIDi[c]];
      const BondDatumSR& bdik = bondDataSR[arik.bond().id()];      
      fdata[0].tM += bdik.M_[1-arik.position()];
    }
    fdata[0].ccArIDs = fullArIDi;
    //build tree
    int k, sz;
    for(c = 0; c < nfracArIDi; c++) {
      const AMod::Arrow& arik = ai.arrows[fracArIDi[c]];
      const BondDatumSR& bdik = bondDataSR[arik.bond().id()];
      sz = fdata.size();
      for(k = 0; k < sz; k++) {
	if(fdata[k].tN+1 < 3) {//new node
	  fdata.push_back(FDatum());
	  fdata.back().W = fdata[k].W*bdik.SdN;
	  fdata.back().tN = fdata[k].tN+1;
	  fdata.back().tM = fdata[k].tM+bdik.M_[1-arik.position()];
	  fdata.back().ccArIDs = fdata[k].ccArIDs;
	  fdata.back().ccArIDs.push_back(fracArIDi[c]);
	}
	fdata[k].W *= (1.0-bdik.SdN);
      }//end for(k...
    }//end for(c...
  }//end if(nfullArIDi < 3)
  if(nfullArIDi+nfracArIDi > 2) {
    double sumW = 0.0;
    int nfdata = fdata.size();
    for(c = 0; c < nfdata; c++) {
      sumW += fdata[c].W;
      fdata[c].tM = min(fdata[c].tM,3.0);
    }
    fdata.push_back(FDatum());
    fdata.back().W = 1.0-sumW;
    fdata.back().tN = 3;
    fdata.back().tM = 0.0;
  }
}

double Energy::fFAT(const AMod::Bond& bondSR) {
  const AMod::Arrow& arij = bondSR.arrow(0);
  const AMod::Arrow& arji = bondSR.arrow(1);
  const AMod::Atom& ai = arij.host();
  const AMod::Atom& aj = arji.host();
  int ck, cl;

  /* prepare FAT */
  std::vector<FDatum> fdataij;
  std::vector<FDatum> fdataji;
  fPreF(arij, fdataij);
  fPreF(arji, fdataji);

  /* Main loop for computing Fconj, A and T */
  double Nelij, Nelji;
  const double Nelmin[4] = {4.0/1, 4.0/(1+1), 4.0/(2+1), 4.0/(3+1)};
  const double Nelmax[4] = {4.0-0, 4.0-1, 4.0-2, 4.0-3};
  double Nconj;
  double aFconj;
  double Fconj;
  double Deltael, a, A;
  int k1, k2, l1, l2;
  double tijk[3], tjil[3];
  double T;
  int nfdataij = fdataij.size();
  int nfdataji = fdataji.size();
  Fconj = 0.0;
  A = 0.0;
  T = 0.0;
  for(ck = 0; ck < nfdataij; ck++) {
    for(cl = 0; cl < nfdataji; cl++) {
      const FDatum& fdij = fdataij[ck];
      const FDatum& fdji = fdataji[cl]; 
      Nelij = (4.0-fdij.tM)/(fdij.tN+1.0-fdij.tM);
      Nelji = (4.0-fdji.tM)/(fdji.tN+1.0-fdji.tM);
      if((fdij.tN == 0 && fdji.tN == 0) ||
	 (fdij.tN == 0 && fdji.tN == 3) ||
	 (fdij.tN == 3 && fdji.tN == 0) ||
	 (fdij.tN == 3 && fdji.tN == 3)) {
	Nconj = 0.0;
      }
      else {
	Nconj = 
	  (Nelij+Nelji-Nelmin[fdij.tN]-Nelmin[fdji.tN])/
	  (Nelmax[fdij.tN]+Nelmax[fdji.tN]-Nelmin[fdij.tN]-Nelmin[fdji.tN]);
      }
      /* Fconj */
      aFconj = 
	(1.0-Nconj)*Const::Fconj0[fdij.tN][fdji.tN]+
	(    Nconj)*Const::Fconj1[fdij.tN][fdji.tN];
      Fconj += fdij.W*fdji.W*aFconj;
      /* A */
      if((fdij.tN == 1 && fdji.tN == 1) ||
	 (fdij.tN == 1 && fdji.tN == 2) ||
	 (fdij.tN == 2 && fdji.tN == 1) ||
	 (fdij.tN == 2 && fdji.tN == 2)) {
	Deltael = Nelij-Nelji;
	a = Const::alpha0*Deltael*Deltael/(1.0+10.0*std::fabs(Deltael));
	A += fdij.W*fdji.W*a;
      }
      /* T */
      if(fdij.tN == 2 && fdji.tN == 2) {
	/* tijk & tjil */
	k1 = fdij.ccArIDs[0];
	k2 = fdij.ccArIDs[1];
	ftvec(arij.dx,ai.arrows[k1].dx,ai.arrows[k2].dx,tijk);
	l1 = fdji.ccArIDs[0];
	l2 = fdji.ccArIDs[1];
	ftvec(arji.dx,aj.arrows[l1].dx,aj.arrows[l2].dx,tjil);
	/* t */
	Deltael = Nelij-Nelji;
	T += fdij.W*fdji.W*ft(costh(tijk,tjil),Nconj,Deltael);
      }
    }//end for(cl...
  }//end for(ck...
  return Fconj+A+T;
}

namespace {

double ft(double y, double z, double Delta) {
  double y2 = y*y;
  double temp = square(z-0.125);
  if(z <= 0.125) {
    double tau1 = Const::At*temp;
    return tau1*square(y2*(1.0-y2));
  }
  double Delta2 = Delta*Delta;
  double tau2 = 
    Const::Bt1*temp*
    square(z+Const::Bt2*Delta2*(Delta2-4.0/9.0))*
    (1.0-Const::Bt3*z)/(Const::Bt4+temp);
  return tau2*(1.0-y2)*square(2.0-y2);
}

void ftvec(const double rij[3], 
	   const double rik1[3], 
	   const double rik2[3], 
	   double t[3]) {
  double wpijk[3], wmijk[3], hwmijk[3];
  double xp[3], xm[3];
  double hrij[3], hrik1[3], hrik2[3];
  double dm;
  normalize(rij,hrij);
  normalize(rik1,hrik1);
  normalize(rik2,hrik2);
  for(int c = 0; c < 3; c++) {
    wpijk[c] = hrik1[c]+hrik2[c];
    wmijk[c] = hrik1[c]-hrik2[c];
  }
  cross(hrij,wpijk,xp);
  cross(hrij,wmijk,xm);
  normalize(wmijk,hwmijk);
  dm = 0.0;
  for(int c = 0; c < 3; c++) dm += hrij[c]*hwmijk[c];
  dm = Const::SQRT3*dm;
  for(int c = 0; c < 3; c++) t[c] = xm[c]+dm*xp[c];
} 

}/* unnamed */

}/* LCBOPII */
