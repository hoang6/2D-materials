#include "Energy.h"
#include "constants.h"
#include "globalConstants.h"
#include "theta2.h"
#include "radic.h"
#include "pij.h"
#include "g.h"
#include <cmath>
#include <limits>

namespace REBO {

Energy::Energy() { 
  tbtools::constantsGen(); 
}

Energy::Energy(const Energy& energy) { *this = energy; }

Energy::~Energy() {}

Energy& Energy::operator= (const Energy& energy) {
  if(this != &energy) {
    bondData = energy.bondData;
  }
  return *this;
}

void Energy::reset() {
  bondData.clear();
}

AMod::BondRange Energy::bondRange() {
  AMod::BondRange br;
  br.insert(6, 6, 0.0, 2.0);
  br.insert(1, 6, 0.0, 1.8);
  br.insert(1, 1, 0.0, 1.7);
  return br;
}

double Energy::cellLength() { return 2.05; }

void Energy::compute_r0(const AMod::Bond& bond) {
  BondDatum& bd = bondData[bond.id()];
  bd.type[0] = (bond.arrow(0).host()->type == 6 ? 1 : 2);
  bd.type[1] = (bond.arrow(1).host()->type == 6 ? 1 : 2);
  bd.fc = tbtools::fc(bond.r, bd.type[0], bd.type[1]);
}

void Energy::compute_r1(const AMod::Bond& bond) {
  int j, k, n;
  int type1, type2, type3;
  double XDB = 0.0, REG = 0.0, costh, W, G, S, Pij;
  BondDatum& bd = bondData[bond.id()];
  for(k = 0; k < 2; k++) {
    const AMod::Arrow& ar = bond.arrow(k);
    const AMod::Atom& a = ar.host();
    n = a.arrows.size();
    type1 = bd.type[k];
    type2 = bd.type[1-k];
    //---------- N ----------
    bd.N[k][0] = bd.N[k][1] = bd.N[k][2] = 0;
    for(j = 0; j < n; j++) {
      if(j == ar.id()) continue;
      const AMod::Arrow& _ar = a.arrows[j];
      BondDatum& _bd = bondData[_ar.bond().id()];
      type3 = _bd.type[1-_ar.position()];
      //----- N -----
      bd.N[k][type3] += _bd.fc;
    }//end for(j...
    bd.N[k][0] = bd.N[k][1]+bd.N[k][2];
    //---------- F ----------
    bd.F[k] = tbtools::F(bd.N[k][0]);
    //---------- XDB, REG, W, costh, G, S----------
    S = 0.0;
    for(j = 0; j < n; j++) {
      if(j == ar.id()) continue;
      const AMod::Arrow& _ar = a.arrows[j];
      BondDatum& _bd = bondData[_ar.bond().id()];
      type3 = _bd.type[1-_ar.position()];
      //----- XDB -----
      XDB = tbtools::XDB(type1, type2, type3);
      //----- REG -----
      REG = tbtools::REG(type1, type2, type3);
      //----- W -----
      W = REG*exp(XDB*(bond.r-_ar.bond().r));
      //----- costh -----
      costh = (ar.dx[0]*_ar.dx[0]+ar.dx[1]*_ar.dx[1]+ar.dx[2]*_ar.dx[2])/(bond.r*_ar.bond().r);
      //----- G -----
      G = tbtools::G(type1, costh, bd.N[k][0]);
      //----- S -----
      S += G*W*bondData[_ar.bond().id()].fc;
    }
    //---------- Pij ----------
    if(type1 == 2) Pij = 0.0;
    else Pij = tbtools::Pij(type2, bd.N[k][2], bd.N[k][1]);
    //---------- bsigmapi ----------
    bd.bsigmapi[k] = 1.0/sqrt(1.0+S+Pij);
  }//end for(k...
}

void Energy::compute_r2(const AMod::Bond& bond) {  
  int i, j, k, n1, n2;
  int type1, type2;
  double bpi;
  const AMod::Arrow& ar1 = bond.arrow(0);
  const AMod::Arrow& ar2 = bond.arrow(1);
  const AMod::Atom& a1 = ar1.host();
  const AMod::Atom& a2 = ar2.host();
  BondDatum& bd = bondData[bond.id()];
  type1 = bd.type[0];
  type2 = bd.type[1];
  //---------- consts that depend only on type ----------
  bd.Q = tbtools::Q(type1, type2);
  bd.A = tbtools::A(type1, type2);
  bd.alpha = tbtools::alpha(type1, type2);
  for(k = 0; k < 3; k++) {
    bd.B[k] = tbtools::B(type1, type2, k+1);
    bd.beta[k] = tbtools::beta(type1, type2, k+1);
  }//end for(n = 0; n < 3...
  //---------- Nconj ----------
  double sum1 = 0.0;
  n1 = a1.arrows.size();
  for(j = 0; j < n1; j++) {
    if(j == ar1.id()) continue;
    const AMod::Arrow& _ar = a1.arrows[j];
    BondDatum& _bd = bondData[_ar.bond().id()];
    if(_bd.type[1-_ar.position()] == 1) { 
      sum1 += _bd.fc*_bd.F[1-_ar.position()];
    }
  }//end for(j...
  double sum2 = 0.0;
  n2 = a2.arrows.size();
  for(j = 0; j < n2; j++) {
    if(j == ar2.id()) continue;
    const AMod::Arrow& _ar = a2.arrows[j];
    BondDatum& _bd = bondData[_ar.bond().id()];
    if(_bd.type[1-_ar.position()] == 1) {
      sum2 += _bd.fc*_bd.F[1-_ar.position()];
    }
  }//end for(j...
  bd.sumFfc[ar1.position()] = sum1;
  bd.sumFfc[ar2.position()] = sum2;
  bd.Nconj = 1.0+sum1*sum1+sum2*sum2;
  //---------- Tij ----------
  bd.Tij = tbtools::Tij(bd.N[ar1.position()][0], bd.N[ar2.position()][0], bd.Nconj);
  //---------- Fij ----------
  bd.Fij = tbtools::Fij(type1, type2, bd.N[ar1.position()][0], bd.N[ar2.position()][0], bd.Nconj);
  //---------- bDH ----------
  double tot;
  if((type1 != 1 || type2 != 1) || fabs(bd.Tij) < tbtools::EPS_DOUBLE) {
    bd.sum_theta2 = 0.0;
    bd.bDH = 0.0;
  }
  else {
    tot = 0.0;
    for(i = 0; i < n1; i++) {
      if(i == ar1.id()) continue;
      const AMod::Arrow& _ar1 = a1.arrows[i];
      for(j = 0; j < n2; j++) {
	if(j == ar2.id()) continue;
	const AMod::Arrow& _ar2 = a2.arrows[j];
	if(_ar1.bond().id() != _ar2.bond().id()) {
	  tot += 
	    (1-tbtools::theta2(_ar1.dx,_ar2.dx,ar1.dx))*
	    bondData[_ar1.bond().id()].fc*
	    bondData[_ar2.bond().id()].fc;
	}//end if(_ar.pbond...
      }//end(j...
    }//end(i...
    bd.sum_theta2 = tot;
    bd.bDH = bd.Tij*tot;
  }//end else
  //---------- bpi ----------
  bpi = bd.Fij+bd.bDH;
  //---------- bbar ----------
  bd.bbar = bd.bsigmapi[0]+bd.bsigmapi[1]+bpi;
  //---------- VR ----------
  bd.VR = bd.fc*(1+bd.Q/bond.r)*bd.A*exp(-bd.alpha*bond.r);
  //---------- VA ----------
  bd.VApre = 0.0;
  for(i = 0; i < 3; i++) bd.VApre += bd.B[i]*exp(-bd.beta[i]*bond.r);
  bd.VA = (bd.VApre*bd.fc/2.0);
  //---------- E ----------
  bd.E = bd.VR-bd.bbar*bd.VA;
}

double Energy::computeE(const AMod::MolTopo& topo) { 
  int k;
  int nbonds = topo.nbonds();
  double tbPot = 0.0;

  bondData.resize(nbonds);

  //---------- must be done in correct order ----------
  //#pragma omp parallel private(k) reduction(+:tbPot) 
  {
    //#pragma omp for
    for(k = 0; k < nbonds; k++) compute_r0(topo.bond(k));
    //#pragma omp for
    for(k = 0; k < nbonds; k++) compute_r1(topo.bond(k));
    //#pragma omp for
    for(k = 0; k < nbonds; k++) compute_r2(topo.bond(k)); 
    //#pragma omp for
    for(k = 0; k < nbonds; k++) tbPot += bondData[k].E;
  }
  
  return tbPot;
}

void Energy::dumpBondData(std::ostream& fdebug, const AMod::Bond& bond) const {
  int k = bond.id();
  const BondDatum& bd = bondData[k];
  int atomIDs[2] = {bond.arrow(0).host().id(), bond.arrow(1).host().id()};
  fdebug << "Bond " << k << std::endl;
  fdebug << "Atom " << atomIDs[0] << "\t" << atomIDs[1] << std::endl;
  fdebug << "type = " << bd.type[0] << "\t" << bd.type[1] << std::endl;
  fdebug << "N[0] = " << bd.N[0][0] << "\t" << bd.N[0][1] << "\t" << bd.N[0][2] << std::endl;
  fdebug << "N[1] = " << bd.N[1][0] << "\t" << bd.N[1][1] << "\t" << bd.N[1][2] << std::endl;
  fdebug << "Nconj = " << bd.Nconj << std::endl;
  fdebug << "F = " << bd.F[0] << "\t" << bd.F[1] << std::endl;
  fdebug << "Fij = " << bd.Fij << std::endl;
  /* debug */
  int type1 = bd.type[0];
  int type2 = bd.type[1];
  const AMod::Arrow& ar1 = bond.arrow(0);
  const AMod::Arrow& ar2 = bond.arrow(1);
  fdebug << "[Debug Fij]" << std::endl;
  fdebug << "Args = "
	 << type1 << " "
	 << type2 << " "
	 << bd.N[ar1.position()][0] << " "
	 << bd.N[ar2.position()][0] << " "
	 << bd.Nconj << std::endl;
  fdebug << "Value = " 
	 << tbtools::Fij(type1, type2, bd.N[ar1.position()][0], bd.N[ar2.position()][0], bd.Nconj) << std::endl;
  /* gubed */
  fdebug << "Tij = " << bd.Tij << std::endl;
  fdebug << "dx = " << bond.arrow(0).dx[0] << "\t" << bond.arrow(0).dx[1] << "\t" << bond.arrow(0).dx[2] << std::endl;
  fdebug << "r = " << bond.r << std::endl;
  fdebug << "fc = " << bd.fc << std::endl;
  fdebug << "bDH = " << bd.bDH << std::endl;
  fdebug << "bsigmapi = " << bd.bsigmapi[0] << "\t" << bd.bsigmapi[1] << std::endl;
  fdebug << "bbar = " << bd.bbar << std::endl;
  fdebug << "VR = " << bd.VR << std::endl;
  fdebug << "VA = " << bd.VA << std::endl;    
  fdebug << std::endl;
}

} /* REBO */
