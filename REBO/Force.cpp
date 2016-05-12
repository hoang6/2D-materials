#include "Force.h"
#include "constants.h"
#include "theta2.h"
#include "radic.h"
#include "pij.h"
#include "g.h"
#include "../Util/util.h"
#include <cmath>

namespace REBO {

Force::Force(Energy& energy):bondData(energy.bondData) {}

Force::~Force() {}

void Force::compute_f0(const AMod::Arrow& ar0, double* f) {
  //derivative of (*ar0.pbond)~E as a function of (ar0.dx)
  int i, j, k, n;
  int type1, type2, type3;
  const AMod::Bond& bond = ar0.bond();
  const BondDatum& bd = bondData[bond.id()]; 
  int sign[2];
  sign[ar0.position()] = +1;
  sign[1-ar0.position()] = -1;
  
  //------------------------- arrow -------------------------
  double XDB, REG, W, costh, G;
  double dW_dr[3], dcosth_dr[3], dG_dcosth, dS_dr[3], bsigmapi3, dbsigmapi_dr[2][3];
  for(k = 0; k < 2; k++) {
    const AMod::Arrow& ar = bond.arrow(k);
    const AMod::Atom& atom = ar.host();
    n = atom.arrows.size();
    type1 = bd.type[k];
    type2 = bd.type[1-k];
    //------ dS/dr ------
    for(i = 0; i < 3; i++) dS_dr[i] = 0.0;
    for(j = 0; j < n; j++) {
      if(j == ar.id()) continue;
      const AMod::Arrow& _ar = atom.arrows[j];
      const AMod::Bond& _bond = _ar.bond();
      const BondDatum& _bd = bondData[_bond.id()];
      type3 = _bd.type[1-_ar.position()];
      //----- dW/dr -----
      XDB = tbtools::XDB(type1, type2, type3);
      REG = tbtools::REG(type1, type2, type3);
      W = REG*exp(XDB*(bond.r-_bond.r));
      for(i = 0; i < 3; i++) 
	dW_dr[i] = sign[k]*(ar.dx[i]/bond.r)*XDB*W;
      //----- dcosth/dr -----
      costh = 
	(ar.dx[0]*_ar.dx[0]+ar.dx[1]*_ar.dx[1]+ar.dx[2]*_ar.dx[2])/(bond.r*_bond.r);
      for(i = 0; i < 3; i++)
	dcosth_dr[i] = 
	  sign[k]*(_ar.dx[i]/(bond.r*_bond.r)-ar.dx[i]*costh/(bond.r*bond.r));
      //----- dG/dcosth -----
      G = tbtools::G(type1, costh, bd.N[k][0]);
      dG_dcosth = tbtools::dG_dcosth(type1, costh, bd.N[k][0]);
      //----- dS/dr -----
      for(i = 0; i < 3; i++) 
	dS_dr[i] += 
	  (dG_dcosth*dcosth_dr[i]*W+G*dW_dr[i])*_bd.fc;
    }//end for(j...
    //------ dsigmpi/dr ------
    bsigmapi3 = bd.bsigmapi[k]*bd.bsigmapi[k]*bd.bsigmapi[k];
    for(i = 0; i < 3; i++) dbsigmapi_dr[k][i] = -0.5*bsigmapi3*dS_dr[i];
  }//end for(k...
  
  //------------------------- bond -------------------------
  double dbDH_dr[3], dbbar_dr[3], dfc_dr[3], dVR_dr[3], dVA_dr[3];
  //------ dbDH/dr------
  if(bd.bDH == 0.0) dbDH_dr[0] = dbDH_dr[1] = dbDH_dr[2] = 0.0;
  else {
    const AMod::Arrow& ar1 = bond.arrow(0);
    const AMod::Arrow& ar2 = bond.arrow(1);
    const AMod::Atom& a1 = ar1.host();
    const AMod::Atom& a2 = ar2.host();
    int n1 = a1.arrows.size();
    int n2 = a2.arrows.size();
    double tot[3], result[3];
    tot[0] = tot[1] = tot[2] = 0.0;
    for(i = 0; i < n1; i++) {
      if(i == ar1.id()) continue;
      const AMod::Arrow& _ar1 = a1.arrows[i];
      for(j = 0; j < n2; j++) {
	if(j == ar2.id()) continue;
	const AMod::Arrow& _ar2 = a2.arrows[j];
	if(_ar1.bond().id() != _ar2.bond().id()) {
	  tbtools::dtheta2_dr(_ar1.dx,_ar2.dx,ar0.dx,result);
	  for(k = 0; k < 3; k++) 
	    tot[k] += 
	      -result[k]*
	      bondData[_ar1.bond().id()].fc*
	      bondData[_ar2.bond().id()].fc;
	}//end if(_ar.pbond...
      }//end(j...
    }//end(i...
    for(i = 0; i < 3; i++) dbDH_dr[i] = bd.Tij*tot[i];
  }//end else {
  //------ dbbar/dr ------
  for(i = 0; i < 3; i++) 
    dbbar_dr[i] = dbsigmapi_dr[0][i]+dbsigmapi_dr[1][i]+dbDH_dr[i];
  //------ dfc/dr ------
  for(i = 0; i < 3; i++) 
    dfc_dr[i] = tbtools::dfc_dr(bond.r,bd.type[0],bd.type[1])*(ar0.dx[i]/bond.r);
  //------ dVR/dr ------
  for(i = 0; i < 3; i++) 
    dVR_dr[i] = 
      (dfc_dr[i]*(1+bd.Q/bond.r)-bd.fc*bd.Q*ar0.dx[i]/(bond.r*bond.r*bond.r))*
      bd.A*exp(-bd.alpha*bond.r)+bd.VR*(-bd.alpha*ar0.dx[i]/bond.r);
  //------ dVA/dr ------
  for(i = 0; i < 3; i++) {
    dVA_dr[i] = 0.0;
    for(j = 0; j < 3; j++) 
      dVA_dr[i] += bd.B[j]*exp(-bd.beta[j]*bond.r)*(-bd.beta[j]*ar0.dx[i]/bond.r);
    dVA_dr[i] = 0.5*(dfc_dr[i]*bd.VApre+bd.fc*dVA_dr[i]);
  }
  //------ f = dE/dr ------
  for(i = 0; i < 3; i++) 
    f[i] = -dVR_dr[i]+dbbar_dr[i]*bd.VA+bd.bbar*dVA_dr[i];
}

void Force::compute_f1(const AMod::Arrow& ar0, const AMod::Arrow& ar1, double* f) {
  int i, j, k;
  const AMod::Bond& bond = ar0.bond();
  const AMod::Atom& atom = ar0.host();
  const BondDatum& bd = bondData[bond.id()];
  const BondDatum& bd_1 = bondData[ar1.bond().id()];

  //------------------------- arrow -------------------------
  int n, type1, type2, type3;
  double _fc, _dfc_dr, dfc_dr[3], dN_dr[3], XDB, REG, W, costh, G, dW_dr, dcosth_dr;
  double dG_dcosth, dG_dN, dS_dr[3], dPij_dr[3], bsigmapi3, dbsigmapi_dr[3];
  n = atom.arrows.size();
  type1 = bd.type[ar0.position()];
  type2 = bd.type[1-ar0.position()];
  type3 = bd_1.type[1-ar1.position()];
  //------ dfc/dr ------
  _dfc_dr = tbtools::dfc_dr(ar1.bond().r,type1,type3);
  for(i = 0; i < 3; i++)
    dfc_dr[i] = _dfc_dr*(ar1.dx[i]/ar1.bond().r);
  //------ dN[0,1,2]/dr ------
  for(i = 0; i < 3; i++) dN_dr[i] = dfc_dr[i];
  //------ dPij/dr ------
  if(type1 == 1 && type3 == 1) {
    for(i = 0; i < 3; i++)
      dPij_dr[i] = 
	tbtools::dPij_dNic(type2, bd.N[ar0.position()][2], bd.N[ar0.position()][1])*dN_dr[i];
  }
  else if(type1 == 1 && type3 == 2) {
    for(i = 0; i < 3; i++) 
      dPij_dr[i] = 
	tbtools::dPij_dNih(type2, bd.N[ar0.position()][2], bd.N[ar0.position()][1])*dN_dr[i];
  }
  else   
    dPij_dr[0] = dPij_dr[1] = dPij_dr[2] = 0.0;
  //------ dS/dr ------
  for(i = 0; i < 3; i++) dS_dr[i] = 0.0;
  for(j = 0; j < n; j++) {
    if(j == ar0.id()) continue;
    const AMod::Arrow& _ar = atom.arrows[j];
    const AMod::Bond& _bond = _ar.bond();
    const BondDatum& _bd = bondData[_bond.id()];
    type3 = _bd.type[1-_ar.position()];
    XDB = tbtools::XDB(type1, type2, type3);
    REG = tbtools::REG(type1, type2, type3);
    W = REG*exp(XDB*(bond.r-_bond.r));
    costh = 
      (ar0.dx[0]*_ar.dx[0]+ar0.dx[1]*_ar.dx[1]+ar0.dx[2]*_ar.dx[2])/(bond.r*_bond.r);
    G = tbtools::G(type1, costh, bd.N[ar0.position()][0]);
    dG_dN = tbtools::dG_dNit(type1, costh, bd.N[ar0.position()][0]);
    if(j == ar1.id()) {
      _fc = bd_1.fc;
      dG_dcosth = tbtools::dG_dcosth(type1, costh, bd.N[ar0.position()][0]);
      for(i = 0; i < 3; i++) {
	dW_dr = W*XDB*(-_ar.dx[i]/_bond.r);
	dcosth_dr = (ar0.dx[i]/(bond.r*_bond.r)-_ar.dx[i]*costh/(_bond.r*_bond.r));
	dS_dr[i] += 
	  (dG_dcosth*dcosth_dr+dG_dN*dN_dr[i])*W*_fc+
	  G*(dW_dr*_fc+W*dfc_dr[i]);
      }
    }
    else {
      for(i = 0; i < 3; i++)
	dS_dr[i] += dG_dN*dN_dr[i]*W*_bd.fc;
    }
  }//end for(j...
  //------ dbsigmapi/dr ------
  bsigmapi3 = bd.bsigmapi[ar0.position()]*bd.bsigmapi[ar0.position()]*bd.bsigmapi[ar0.position()];
  for(i = 0; i < 3; i++)
    dbsigmapi_dr[i] = -0.5*bsigmapi3*(dS_dr[i]+dPij_dr[i]);
 
  //------------------------- bond -------------------------
  double dNconj_dr[3], dTij_dN, dTij_dNconj, dFij_dN, dFij_dNconj, dbDH_dr[3], dbbar_dr[3];
  //------ dNconj/dr ------
  for(i = 0; i < 3; i++) dNconj_dr[i] = 0.0;
  if(type3 == 1) {
    for(i = 0; i < 3; i++) 
      dNconj_dr[i] = 2*bd.sumFfc[ar0.position()]*dfc_dr[i]*bd_1.F[1-ar1.position()];
  }
  //------ dTij & dFij ------
  dTij_dN = 
    tbtools::dTij_dXNT1(bd.N[ar0.position()][0], bd.N[1-ar0.position()][0], bd.Nconj);
  dTij_dNconj = 
    tbtools::dTij_dCONJUG(bd.N[ar0.position()][0], bd.N[1-ar0.position()][0], bd.Nconj);
  dFij_dN = 
    tbtools::dFij_dXNT1(type1, type2, bd.N[ar0.position()][0], bd.N[1-ar0.position()][0], bd.Nconj);
  dFij_dNconj = 
    tbtools::dFij_dCONJUG(type1, type2, bd.N[ar0.position()][0], bd.N[1-ar0.position()][0], bd.Nconj);
  if(bd.bDH == 0.0) dbDH_dr[0] = dbDH_dr[1] = dbDH_dr[2] = 0.0;
  else {
    const AMod::Arrow& ar2 = ar0.brother();
    const AMod::Atom& a2 = ar2.host();
    int n2 = a2.arrows.size();
    double tot[3], result[3];
    tot[0] = tot[1] = tot[2] = 0.0;
    for(j = 0; j < n2; j++) {
      if(j == ar2.id()) continue;
      const AMod::Arrow& _ar2 = a2.arrows[j];
      if(ar1.bond().id() != _ar2.bond().id()) {
	tbtools::dtheta2_da(ar1.dx,_ar2.dx,ar0.dx,result);
	for(k = 0; k < 3; k++) 
	  tot[k] += 
	    (-result[k]*bd_1.fc+
	     (1-tbtools::theta2(ar1.dx,_ar2.dx,ar0.dx))*dfc_dr[k])*
	    bondData[_ar2.bond().id()].fc;
      }//end if(_ar.pbond...
    }//end(j...
    for(k = 0; k < 3; k++) 
      dbDH_dr[k] = 
	(dTij_dN*dN_dr[k]+dTij_dNconj*dNconj_dr[k])*bd.sum_theta2+
	bd.Tij*tot[k];
  }//end else {
  for(i = 0; i < 3; i++) 
    dbbar_dr[i] = 
      dbsigmapi_dr[i]+
      dFij_dN*dN_dr[i]+dFij_dNconj*dNconj_dr[i]+dbDH_dr[i];
  for(i = 0; i < 3; i++)
    f[i] = dbbar_dr[i]*bd.VA;
}

void Force::compute_f2(const AMod::Arrow& ar0, const AMod::Arrow& ar1, const AMod::Arrow& ar2, double* f) {
  int i;
  const AMod::Bond& bond = ar0.bond();
  const AMod::Bond& bond_1 = ar1.bond();
  const AMod::Bond& bond_2 = ar2.bond();
  const BondDatum& bd = bondData[bond.id()];
  const BondDatum& bd_1 = bondData[bond_1.id()];
  const BondDatum& bd_2 = bondData[bond_2.id()];

  double _dfc_dr, dfc_dr[3], dN_dr[3], dF_dN, dNconj_dr[3];
  double dTij_dNconj, dFij_dNconj, dbDH_dr[3], dbbar_dr[3];
  _dfc_dr = tbtools::dfc_dr(bond_2.r, bd_2.type[0], bd_2.type[1]);
  for(i = 0; i < 3; i++) 
    dfc_dr[i] = _dfc_dr*ar2.dx[i]/bond_2.r;
  for(i = 0; i < 3; i++) dN_dr[i] = dfc_dr[i];
  dF_dN = tbtools::dF_dN(bd_1.N[1-ar1.position()][0]);
  for(i = 0; i < 3; i++) 
    dNconj_dr[i] = 2*bd.sumFfc[ar0.position()]*bd_1.fc*dF_dN*dN_dr[i];
  dTij_dNconj = 
    tbtools::dTij_dCONJUG(bd.N[ar0.position()][0], bd.N[1-ar0.position()][0], bd.Nconj);
  dFij_dNconj = 
    tbtools::dFij_dCONJUG(bd.type[ar0.position()], bd.type[1-ar0.position()], 
			  bd.N[ar0.position()][0], bd.N[1-ar0.position()][0], bd.Nconj);
  if(bd.bDH == 0.0) dbDH_dr[0] = dbDH_dr[1] = dbDH_dr[2] = 0.0;
  else {
    for(i = 0; i < 3; i++) dbDH_dr[i] = dTij_dNconj*dNconj_dr[i]*bd.sum_theta2;
  }
  for(i = 0; i < 3; i++) 
    dbbar_dr[i] = dFij_dNconj*dNconj_dr[i]+dbDH_dr[i];
  for(i = 0; i < 3; i++)
     f[i] = dbbar_dr[i]*bd.VA;
}

void Force::computeBF(const AMod::Arrow& ar, double* f) {
  int i, imax, j, jmax, k;
  int x, sign;
  double f_temp[3];
  //self
  compute_f0(ar, f);
  //1st and 2nd neighb
  for(k = 0; k < 2; k++) {
    sign = (k == 0 ? +1 : -1);
    const AMod::Arrow& ar0 = (k == 0 ? ar : ar.brother());
    const AMod::Atom& a0 = ar0.host();
    //1st neighb
    for(j = 0, jmax = a0.arrows.size(); j < jmax; j++) {
      if(j == ar0.id()) continue;
      const AMod::Arrow& ar1 = a0.arrows[j];
      compute_f1(ar1, ar0, f_temp);
      for(x = 0; x < 3; x++) f[x] += sign*f_temp[x];
      //2nd neighb
      const AMod::Arrow& ar1b = ar1.brother();
      const AMod::Atom& a1 = ar1b.host();
      for(i = 0, imax = a1.arrows.size(); i < imax; i++) {
	if(i == ar1b.id()) continue;
	const AMod::Arrow& ar2 = a1.arrows[i];
	compute_f2(ar2, ar1b, ar0, f_temp);
	for(x = 0; x < 3; x++) f[x] += sign*f_temp[x];
      }//end for(i...
    }//end for(j...
  }//end for(k...
}

void Force::computeAF(const AMod::Atom& atom, double* f) {
  int x, k, n = atom.arrows.size();
  double f_temp[3];
  f[0] = f[1] = f[2] = 0.0;
  for(k = 0; k < n; k++) {
    computeBF(atom.arrows[k], f_temp);
    for(x = 0; x < 3; x++) f[x] += -f_temp[x];
  }   
}

void Force::computeF(AMod::Molecule& mol) {
  int j, jmax, k, n, x, sign, naxes, natoms, nbonds;
  double* f;
  naxes = mol.naxes();
  natoms = mol.natoms();
  nbonds = mol.nbonds();
  //calculate bond forces
  //#pragma omp parallel private(k)
  {
    //#pragma omp for
    for(k = 0; k < nbonds; k++) {
      computeBF(mol.bond(k).arrow(0), bondData[k].f);
    }
  }
  //init axes force
  for(k = 0; k < naxes; k++)
    Util::vassign(0.0, mol.axis(k).fpot, 3);
  //calculate atom forces from bond forces
  for(k = 0; k < natoms; k++) {
    AMod::Atom& atom = mol.atom(k);
    Util::vassign(0.0, atom->fpot, 3);
    jmax = atom.arrows.size();
    for(j = 0; j < jmax; j++) {
      AMod::Arrow& ar = atom.arrows[j];
      sign = (ar.position() == 0 ? +1 : -1);
      f = bondData[ar.bond().id()].f;
      for(x = 0; x < 3; x++)
	atom->fpot[x] += -sign*f[x];
      for(n = 0; n < naxes; n++) {
	if(ar.n[n] > 0) {
	  for(x = 0; x < 3; x++)
	    mol.axis(n).fpot[x] += sign*f[x]*ar.n[n];
	}
      }//end for(n...
    }//end for(j... 
  }//end for(k... 
}

}/* REBO */
