#include "MolAdjuster.h"
#include "MolAnal.h"
#include "Molecule.h"
#include "../Util/util.h"
#include "../NR/jacobi.h"
#include "../NR/eigsrt.h"

namespace AMod {

MolAdjuster::MolAdjuster() {}

MolAdjuster::MolAdjuster(Molecule& _mol): MolServ(_mol) {}

void MolAdjuster::reset() {}

void MolAdjuster::adjustMassCenter(const Point& _massCenter) const {
  Point curr_mc = massCenter();
  move(curr_mc, _massCenter);
}

MolAdjuster::Point MolAdjuster::massCenter(const Molecule& _mol){
  int j, k, n = _mol.natoms();
  double mass_k, mass_tot = 0.0;
  Point curr_mc;
  for(k = 0; k < n; k++) {
    const Atom& a = _mol.atom(k);
    mass_k = PTE::atomicMass(a->type);
    mass_tot += mass_k;
    for(j = 0; j < 3; j++) curr_mc[j] += a->x[j]*mass_k;
  }//end for(k...
  for(j = 0; j < 3; j++) curr_mc[j] /= mass_tot;
  return curr_mc;
}

void MolAdjuster::adjustPrinpAxes(const Axes& _prinpAxes) const {
  Point curr_mc = massCenter();
  adjustAxes(molecule(), prinpAxes(), _prinpAxes);
  adjustMassCenter(curr_mc);
}

MolAdjuster::Axes MolAdjuster::prinpAxes(const Molecule& _mol){
  NR::Matrix<double> mI(1,3,1,3);
  NR::Matrix<double> mv(1,3,1,3);
  NR::Vector<double> vd(1,3);
  double** I = mI.pointer();
  double** v = mv.pointer();
  double* d = vd.pointer();
  Axes curr_pa;
  int j, k, nrot;
  //prepare inertial matrix
  inertial(_mol, I[1][1],I[2][2],I[3][3],I[1][2],I[2][3],I[1][3]);
  I[2][1] = I[1][2];
  I[3][2] = I[2][3];
  I[1][3] = I[3][1];
  //find eig values and eig vectors
  NR::jacobi(I,3,d,v,&nrot);
  //copy to paxes
  for(j = 0; j < 3; j++)
    for(k = 0; k < 3; k++)
      curr_pa[k][j] = v[j+1][k+1];
  return curr_pa;
}

void MolAdjuster::adjust(const Point& _massCenter, 
			 const Axes& _prinpAxes) const {
  adjustAxes(molecule(), prinpAxes(), _prinpAxes);
  adjustMassCenter(_massCenter);
}

void MolAdjuster::adjust(const Molecule& _mol) const {
  adjust(massCenter(_mol), prinpAxes(_mol));
}

void MolAdjuster::move(const Point& _oldPt, const Point& _newPt) const {
  adjustPosition(molecule(), _oldPt, _newPt);
  molecule().topo().update(CASE_TRANSLATION);
}

void MolAdjuster::rotate(const Axes& _oldAxes, const Axes& _newAxes) const {
  adjustAxes(molecule(), _oldAxes, _newAxes);
  molecule().topo().update(CASE_ROTATION);
}

bool MolAdjuster::align(Molecule& mol, int bondID, 
			const Molecule& _mol, int _bondID) {
  Point curr_mc;
  const Bond& _bond = _mol.bond(_bondID);
  Bond& bond = mol.bond(bondID);
  Axes _axes, axes;
  Util::vcopy(_bond.arrow(0).dx, _axes[0], 3);
  Util::vcopy(bond.arrow(0).dx, axes[0], 3);
  if(!MolAnal::bondNorm(_bond, _axes[2], NULL)) return false;
  if(!MolAnal::bondNorm(bond, axes[2], NULL)) return false;
  Util::cross3(_axes[2],_axes[0],_axes[1]); 
  Util::cross3(_axes[0],_axes[1],_axes[2]);
  Util::cross3(axes[2],axes[0],axes[1]); 
  Util::cross3(axes[0],axes[1],axes[2]);
  for(int k = 0; k < 3; k++) {
    if(Util::vdotv(_axes[k],_axes[k],3) < 1e-16) return false;
    else Util::normalize(_axes[k],3);
    if(Util::vdotv(axes[k],axes[k],3) < 1e-16) return false;
    else Util::normalize(axes[k],3);
  }
  if(Util::vdotv(_axes[2],axes[2],3) < 0) {
    Util::negv(_axes[1],3);
    Util::negv(_axes[2],3);
  }
  //rot 90
  for(int k = 0; k < 3; k++) Util::exch(_axes[0][k],_axes[1][k]);
  Util::negv(_axes[1],3);
  if(Util::vdotv(_axes[0],axes[0],3) < 0) {
    Util::negv(_axes[0],3);
    Util::negv(_axes[1],3);
  }
  //move & rotate
  curr_mc = massCenter(mol);
  adjustAxes(mol, axes, _axes);
  adjustPosition(mol, massCenter(mol), curr_mc);
  mol.topo().update(CASE_ROTATION);
  return true;
}

void MolAdjuster::inertial(const Molecule& _mol,
			   double& Ixx, double& Iyy, double& Izz,
			   double& Ixy, double& Iyz, double& Izx) {
  //prepare inertial tensor I
  int k, n = _mol.natoms();
  double mass_k;
  Ixx = Iyy = Izz = Ixy = Iyz = Izx = 0.0;
  for(k = 0; k < n; k++) {
    const Atom& a = _mol.atom(k);
    mass_k = PTE::atomicMass(a->type);
    Ixx += +mass_k*(a->x[1]*a->x[1]+a->x[2]*a->x[2]);
    Ixy += -mass_k*a->x[0]*a->x[1];
    Izx += -mass_k*a->x[0]*a->x[2];
    Iyy += +mass_k*(a->x[2]*a->x[2]+a->x[0]*a->x[0]);
    Iyz += -mass_k*a->x[1]*a->x[2];
    Izz += +mass_k*(a->x[0]*a->x[0]+a->x[1]*a->x[1]);
  }//end for(k...
}

void MolAdjuster::adjustPosition(Molecule& _mol,
				 const Point& _oldPt, 
				 const Point& _newPt) {
  int j, k, n = _mol.natoms();
  for(k = 0; k < n; k++)
    for(j = 0; j < 3; j++)
      _mol.atom(k)->x[j] += _newPt[j]-_oldPt[j];
}

void MolAdjuster::adjustAxes(Molecule& _mol,
			     const Axes& _oldAxes, 
			     const Axes& _newAxes) {
  int j, k, nax, nt;
  double a[3];
  Atom::Data* pdata;
  nax = _mol.naxes();
  nt = nax+_mol.natoms();
  for(k = 0; k < nt; k++) {
    if(k < nax) 
      pdata = &_mol.axis(k);
    else
      pdata = &_mol.atom(k-nax).data();
    for(j = 0; j < 3; j++) 
      a[j] = Util::vdotv(pdata->x,_oldAxes[j],3);
    for(j = 0; j < 3; j++) 
      pdata->x[j] = a[0]*_newAxes[0][j]+a[1]*_newAxes[1][j]+a[2]*_newAxes[2][j];
  }//end for(k = 0...
}

} /* AMod */
