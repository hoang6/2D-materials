#include "MolFactory.h"
#include "Molecule.h"
#include "C60XYZ.h"
#include "constants.h"
#include "../Util/util.h"
#include "../NR/util.h"
#include "../Util/constants.h"
#include <cmath>

namespace AMod {

MolFactory::MolFactory() {}

MolFactory::MolFactory(Molecule& _mol): MolServ(_mol) {}

void MolFactory::reset() {}

void MolFactory::makeGPH_PBC(const int chiV[2], const int heightV[2]) {
  makeGPH(chiV, heightV);
  molecule().setStdAxes();
  molecule().setFixedAxes(2);
}

void MolFactory::makeGPH_NPBC(const int chiV[2], const int heightV[2], bool _capH) {
  makeGPH_PBC(chiV, heightV);
  clearBoundaryBonds();
  if(_capH) capHydrogen();
  molecule().topo().update(CASE_AXES);
  molecule().setFixedAxes(0);
}

void MolFactory::makeCNT_PBC(const int chiV[2], const int heightV[2]) {
  makeCNT(chiV, heightV);
  molecule().setStdAxes();
  molecule().setFixedAxes(1);
}

void MolFactory::makeCNT_NPBC(const int chiV[2], const int heightV[2], bool _capH) {
  makeCNT_PBC(chiV, heightV);
  clearBoundaryBonds();
  if(_capH) capHydrogen();
  molecule().topo().update(CASE_AXES);
  molecule().setFixedAxes(0);
}

void MolFactory::makeC60() {
  Molecule& mol = molecule();
  mol.axes().clear();
  mol.atomsData().resize(60);
  for(int k = 0; k < 60; k++) {
    for(int j = 0; j < 3; j++) {
      mol.atomsData()[k].type = PTE::Carbon;
      mol.atomsData()[k].x[j] = C60XYZ[3*k+j];
    }
  }
  mol.topo().sync();
  mol.topo().update(CASE_AXES);
  mol.topo().update(CASE_BRUTAL);
  molecule().setStdAxes();
  molecule().setFixedAxes(0);
}

void MolFactory::makeGPH(const int chiV[2], const int heightV[2]) {
  Molecule& mol = molecule();
  int col, row, k;
  int colmin, rowmin;
  int colmax, rowmax;
  int count;
  int id,nb_c,nb_r;
  int** aidMat[2];
  double x,y;
  Box box;
  const double a = Box::SCALE;
  const double dr[2][2] = {{0.5, +0.25/Util::SIN60},
			   {0.5, -0.25/Util::SIN60}};

  //--------------- initialize ---------------
  box.prepare(chiV,heightV);
  boxInfo = box.information();
  colmin = boxInfo.colRange[0];
  rowmin = boxInfo.rowRange[0];
  colmax = boxInfo.colRange[1];
  rowmax = boxInfo.rowRange[1];
  aidMat[0] = NR::matrix<int>(colmin,colmax,rowmin,rowmax);
  aidMat[1] = NR::matrix<int>(colmin,colmax,rowmin,rowmax);
  
  //--------------- create atoms ---------------
  /*** atoms ***/
  mol.axes().resize(2);
  mol.atomsData().resize(boxInfo.nGrids*2);
  count = 0;
  for(col = colmin; col <= colmax; col++) {
    for(row = rowmin; row <= rowmax; row++) {
      if(box.inQuad(col,row)) {
	box.cr2xy(col,row,x,y);
	for(k = 0; k < 2; k++) {
	  Atom::Data& data = mol.atomsData()[count++];
	  data.type = PTE::Carbon;
	  data.x[0] = a*(x+dr[k][0]);
	  data.x[1] = a*(y+dr[k][1]);
	  data.x[2] = 0.0;
	  Util::rot(-boxInfo.chi, data.x[0], data.x[1]);
	  aidMat[k][col][row] = data.id();
	}	
      }//end if(box.inQuad...
    }//end for(row...
  }//end for(col...
  /*** axes ***/
  for(k = 0; k < 2; k++) {
    mol.axis(k).x[0] = a*boxInfo.dQuad[k][0];
    mol.axis(k).x[1] = a*boxInfo.dQuad[k][1];
    mol.axis(k).x[2] = 0.0;
    Util::rot(-boxInfo.chi, mol.axis(k).x[0], mol.axis(k).x[1]);
  }
  /*** topo ***/
  mol.topo().sync();
  mol.topo().update(CASE_AXES);
  mol.topo().update(CASE_LONELY);

  //--------------- create bonds---------------
  for(col = colmin; col <= colmax; col++) {
    for(row = rowmin; row <= rowmax; row++) {
      if(box.inQuad(col,row)) {
	id = aidMat[1][col][row];
	/*** 1st bond ***/
	mol.topo().newBond(id, aidMat[0][col][row]);
	/*** 2nd bond ***/
	box.neighbor(col,row,Box::BOTTOM_LEFT,nb_c,nb_r);
	box.move2quad(nb_c,nb_r,nb_c,nb_r);
	mol.topo().newBond(id, aidMat[0][nb_c][nb_r]);
	/*** 3rd bond ***/
	box.neighbor(col,row,Box::BOTTOM_RIGHT,nb_c,nb_r);
	box.move2quad(nb_c,nb_r,nb_c,nb_r);
	mol.topo().newBond(id, aidMat[0][nb_c][nb_r]);
      }//end if(box.inQuad...
    }//end for(row...
  }//end for(col...

  //--------------- finalize ---------------
  NR::free_matrix<int>(aidMat[0],colmin,colmax,rowmin,rowmax);
  NR::free_matrix<int>(aidMat[1],colmin,colmax,rowmin,rowmax);
}

void MolFactory::makeCNT(const int chiV[2], const int heightV[2]) {
  Molecule& mol = molecule();
  double r, lz;

  makeGPH(chiV,heightV);
  r = mol.axis(0).x[0]/2/Util::PI;
  lz = mol.axis(1).x[1];

  //--------------- axes ---------------
  mol.axes().resize(1);
  mol.axis(0).reset();
  mol.axis(0).x[2] = lz;
  //--------------- atoms ---------------
  for(int k = 0, kmax = mol.natoms(); k < kmax; k++) {
    rollup(r, mol.atom(k)->x);
  }
  //--------------- topo ---------------
  mol.topo().update(CASE_ROTATION);
}

void MolFactory::clearBoundaryBonds() const {
  Molecule& mol = molecule();
  int* n;
  for(int k = 0; k < mol.nbonds(); k++) {
    n = mol.bond(k).arrow(0).n;
    if((n[0]!=0 || n[1]!=0 || n[2]!=0)) {
      mol.topo().deleteBond(k--);
    }
  }//end for(k...
}

void MolFactory::capHydrogen() const {
  Molecule& mol = molecule();
  MolTopoMode mode = mol.topo().mode();
  int j, k, natoms = mol.natoms();
  double tx[3];
  Atom::Data data;
  data.type = PTE::Hydrogen;
  mol.topo().mode(MODE_FREEZE);
  mol.topo().update(CASE_CHANGE_MODE);
  for(k = 0; k < natoms; k++) {
    Atom& a = mol.atom(k);
    if(a->type != PTE::Carbon) continue;
    if(a.arrows.size() == 1) {//change C to H 
      Arrow& ar = a.arrows[0].brother();
      for(j = 0; j < 3; j++)
	tx[j] = ar.host()->x[j]+ar.dx[j]*CH_BOND_LENGTH/CC_BOND_LENGTH;
      a->type = PTE::Hydrogen;
      mol.moveAtom(k,tx);
    }//end if(a.nbonds() == 1)
    else if(a.arrows.size() == 2) {//cap H on C
      for(j = 0; j < 3; j++)
	tx[j] = a.arrows[0].dx[j]+a.arrows[1].dx[j];
      Util::kdotv(CH_BOND_LENGTH/Util::norm2(tx,3),tx,tx,3);
      for(j = 0; j < 3; j++) data.x[j] = a->x[j]-tx[j];
      mol.topo().newBond(a.id(), mol.newAtom(-1,data).id());
    }//end else if(a.nbonds() == 2)
  }//end for(k...
  mol.topo().mode(mode);
  mol.topo().update(CASE_CHANGE_MODE);
}

/*
** Find (x,y)-|(c,r) in coord a=(d,0), b=(d/2,sqrt(3)*d/2) where d=C-C*sqrt(3)
** i)  x/y = k = -(c+2*r)/(2*c+r) where x,y,c,r are integers
** ii) h^2 = d^2*(x^2+y^2+x*y);
** If x = a*tx and y = a*ty where a, tx, ty are integers
** then h = d*a*sqrt(tx^2+ty^2+tx*ty) ==> a = ...
 */
double MolFactory::findHeight(int chic, int chir, double height) {
  int tx = +(chic+2*chir);
  int ty = +(2*chic+chir);
  int k, kmax;
  double f, a;
  //find gcd of tx & ty
  kmax = (tx<ty ? tx : ty);
  for(k = kmax; k >= 1; k--) {
    if(tx%k==0 && ty%k==0) break;
  }
  tx = -(tx/k);
  ty = +(ty/k);
  f = std::sqrt(double(tx*tx+ty*ty+tx*ty));
  a = std::floor(height/f+0.5);
  return f*a;
}

void MolFactory::rollup(const double& r, double* xyz) {
  double th = xyz[0]/r;
  double z = xyz[1];
  xyz[0] = r*cos(th);
  xyz[1] = r*sin(th);
  xyz[2] = z;
}

} /* AMod */
