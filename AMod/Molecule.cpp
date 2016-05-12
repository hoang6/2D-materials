#include "Molecule.h"
#include "MolAdjuster.h"
#include "../Util/util.h"
#include <cmath>

namespace AMod {

Molecule::Molecule() { 
  attach();
  reset();
}

Molecule::Molecule(const Molecule& mol) {
  attach();
  *this = mol;
}

Molecule& Molecule::operator= (const Molecule& mol) {
  if(this != &mol) {
    axes = mol.axes;
    atomsData = mol.atomsData;
    potential = mol.potential;
    work = mol.work;
    io = mol.io; 
    io().attach(*this);
    topo().sync();
    topo = mol.topo;
  }
  return *this;
}

void Molecule::reset() {
  axes().clear();
  atomsData().clear();
  potential(0.0);
  work(0.0);
  topo().reset();
  io().reset();
}

void Molecule::clear() {
  axes().clear();
  atomsData().clear();
  potential(0.0);
  work(0.0);
  topo().clear();
}

Atom& Molecule::newAtom(int _atomID, const Atom::Data& _atomData) {  
  if(_atomID < 0)
    return topo().newAtom(_atomID, atomsData().insert(_atomData));
  return topo().newAtom(_atomID, atomsData().insert(_atomID,_atomData));
}

void Molecule::deleteAtom(int _atomID) {
  atomsData().erase(_atomID);
  topo().deleteAtom(_atomID);
}
 
Atom& Molecule::moveAtom(int _atomID, const double* _x) {
  double* x = atomsData()[_atomID].x;
  for(int k = 0; k < 3; k++) x[k] = _x[k];
  return topo().moveAtom(_atomID);
}

Atom& Molecule::newLonelyAtom(int _atomID, const AData& _atomData) {
  if(_atomID < 0)
    return topo().newLonelyAtom(_atomID, atomsData().insert(_atomData));
  return topo().newLonelyAtom(_atomID, atomsData().insert(_atomID,_atomData));
}

Atom& Molecule::moveLonelyAtom(int _atomID, const double* _x) {
  double* x = atomsData()[_atomID].x;
  for(int k = 0; k < 3; k++) x[k] = _x[k];
  return topo().moveLonelyAtom(_atomID);
}

Atom& Molecule::makeLonelyAtom(int _atomID) {
  return topo().makeLonelyAtom(_atomID);
}

void Molecule::setStdAxes() {
  const double delta = 20.0;
  const double scale = 10.0;
  int k, n = naxes();
  double rmin = 0.0;
  double rmax = 0.0;
  double* tempAx[3];
  MolAdjuster::Axes oldAx, newAx;

  axes().resize(3);

  if(n >= 3) return;

  //add new box edges
  for(k = 0; k < 3; k++) tempAx[k] = axis(k).x;
  Util::build_axes(n, tempAx);
  for(k = n; k < 3; k++) {
    getRange(axis(k).x, rmin, rmax);
    Util::kdotv((rmax-rmin+delta)*scale, axis(k).x, axis(k).x, 3);
    Util::vassign(1,axis(k).fixed,3);
  }
  if(n == 1) {
    axes().exch(0,1);
    axes().exch(1,2);
  }
  topo().update(CASE_AXES);

  //adjust box position
  newAx = MolAdjuster::Axes(1.0,0.0);
  if(n == 2) {
    Util::vcopy(axis(2).x,oldAx[2],3);
    Util::vcopy(axis(0).x,oldAx[0],3);
    Util::normalize(oldAx[2],3);
    Util::normalize(oldAx[0],3);
    Util::cross3(oldAx[2],oldAx[0],oldAx[1]);
    MolAdjuster(*this).rotate(oldAx,newAx);
  }
  else if(n == 1) {
    Util::vcopy(axis(0).x,oldAx[0],3);
    Util::vcopy(axis(1).x,oldAx[1],3);
    Util::vcopy(axis(2).x,oldAx[2],3);
    Util::normalize(oldAx[0],3);
    Util::normalize(oldAx[1],3);
    Util::normalize(oldAx[2],3);
    MolAdjuster(*this).rotate(oldAx,newAx);
  }
}

void Molecule::setFixedAxes(int _dim) {
  if(_dim == 3) {
    for(int k = 0; k < 3; k++)
      Util::vassign(0, axis(k).fixed, 3);
  }
  else if(_dim == 2) {
    axis(0).fixed[0] = 0; axis(0).fixed[1] = 1; axis(0).fixed[2] = 1;
    axis(1).fixed[0] = 0; axis(1).fixed[1] = 0; axis(1).fixed[2] = 1;
    axis(2).fixed[0] = 1; axis(2).fixed[1] = 1; axis(2).fixed[2] = 1;
  }
  else if(_dim == 1) {
    axis(0).fixed[0] = 1; axis(0).fixed[1] = 1; axis(0).fixed[2] = 1;
    axis(1).fixed[0] = 1; axis(1).fixed[1] = 1; axis(1).fixed[2] = 1;
    axis(2).fixed[0] = 1; axis(2).fixed[1] = 1; axis(2).fixed[2] = 0;
  }
  else if(_dim == 0) {
    for(int k = 0; k < 3; k++)
      Util::vassign(1, axis(k).fixed, 3);
  }
}

int Molecule::countAtoms(int _atomType) const {
  int count = 0;
  for(int k = 0, kmax = natoms(); k < kmax; k++) 
    if(atom(k)->type == _atomType) count++;
  return count;
}

void Molecule::axes2array(double* _h) const {
  int k, kmax = naxes();
  for(k = 0; k < kmax; k++) 
    Util::vcopy(axis(k).x, _h+3*k, 3);
}

double Molecule::volume() const {
  return Util::abs(Util::cross_dot3(axis(0).x,axis(1).x,axis(2).x));
}

void Molecule::attach() {
  io().attach(*this);
  topo().attach(*this);
}

void Molecule::getRange(const double _n[], double& _rmin, double& _rmax) const {
  bool flag = true;
  double val;
  for(int k = 0, kmax = natoms(); k < kmax; k++) {
    val = Util::vdotv(_n, atom(k)->x, 3);
    if(flag) { _rmin = val; _rmax = val; flag = false; }
    else {
      _rmin = Util::min(_rmin, val);
      _rmax = Util::max(_rmax, val);
    }
  }
}

} /* AMod */
