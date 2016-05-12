#include "MCAgent_STube.h"
#include <cmath>

namespace MC {

MCAgent_STube::MCAgent_STube() {}

MCAgent_STube::MCAgent_STube(AMod::Molecule& _mol): MCAgent(_mol) {}

MCAgent_STube::~MCAgent_STube() {}

void MCAgent_STube::init() {
  AMod::Molecule& mol = molecule();
  int i, natoms = mol.natoms();
  double score = 0.0, tscore = 0.0;
  _r = 0.0;
  _lz = mol.axis(2).x[2];
  _th.resize(natoms);
  for(i = 0; i < natoms; i++) {
    const AMod::Atom& atom = mol.atom(i);
    _r += std::sqrt(Util::square(atom->x[0])+Util::square(atom->x[1]));
    _th[i] = std::atan2(atom->x[1],atom->x[0]);
    tscore = Util::abs(atom->x[2]-_lz/2.0);
    if(i == 0 || score > tscore) {
      _thFixID = i;
      score = tscore;
    }
  }
  _r = _r/natoms;  
  for(i = 0; i < natoms; i++) {
    AMod::Atom& atom = mol.atom(i);
    atom->x[0] = _r*cos(_th[i]);
    atom->x[1] = _r*sin(_th[i]);
  }
}

void MCAgent_STube::pertAtom(double _dxmax, int _id, double _new_x[3]) {
  const AMod::Atom& atom = molecule().atom(_id);
  _new_id = _id;
  _new_th = _th[_id];
  if(_id != _thFixID)
    _new_th += (_dxmax/_r)*urand(-1,1);
  _new_x[0] = _r*cos(_new_th);
  _new_x[1] = _r*sin(_new_th);
  _new_x[2] = atom->x[2]+_dxmax*urand(-1,1);
}

void MCAgent_STube::accPertAtom(bool _flag) {
  if(_flag)
    _th[_new_id] = _new_th;
}

void MCAgent_STube::pertAxes(double _dhmax) {
  AMod::Molecule& mol = molecule();
  int i, natoms = mol.natoms();

  //old_r & old_lz

  //new_r & new_lz
  _new_r = _r+_dhmax*urand(-1,1);
  _new_lz = _lz+_dhmax*urand(-1,1);
  
  //perturb axes
  mol.axis(2).x[2] = _new_lz;
  
  //perturb atoms
  for(i = 0; i < natoms; i++) {
    mol.atom(i)->x[0] *= _new_r/_r;
    mol.atom(i)->x[1] *= _new_r/_r;
    mol.atom(i)->x[2] *= _new_lz/_lz;
  }
}

void MCAgent_STube::accPertAxes(bool _flag) {
  if(_flag) {
    _r = _new_r;
    _lz = _new_lz;
  }
}

}/* MC */
