#include "MCAgent.h"
#include "../AMod/Molecule.h"

namespace MC {

MCAgent::MCAgent() {}

MCAgent::MCAgent(AMod::Molecule& _mol): AMod::MolServ(_mol) {}

MCAgent::~MCAgent() {}

void MCAgent::init() {}

void MCAgent::final() {}

void MCAgent::pertAtom(double _dxmax, int _id, double _new_x[3]) {
  const AMod::Atom& atom = molecule().atom(_id);
  for(int k = 0; k < 3; k++) {
    if(atom->fixed[k]) 
      _new_x[k] = atom->x[k];
    else
      _new_x[k] = atom->x[k]+_dxmax*urand(-1,1);
  }
}

void MCAgent::accPertAtom(bool _flag) {}

void MCAgent::pertAxes(double _dhmax) {
  AMod::Molecule& mol = molecule();
  int i, j, natoms = mol.natoms();
  double d_h[9];
  double old_h[9], inv_old_h[9];
  double new_h[9], inv_new_h[9];
  double s_x[3];

  //old_h & inv_old_h
  mol.axes2array(old_h);
  Util::minverse(old_h, inv_old_h, 3);
 
  //dh
  Util::vassign(0.0,d_h,9);
  for(i = 0; i < 3; i++) {
    for(j = i; j < 3; j++) {
      if(!mol.axis(i).fixed[j] && !mol.axis(j).fixed[i])
	d_h[3*i+j] = _dhmax*urand(-1,1);
      if(i != j) d_h[3*j+i] = d_h[3*i+j];
    }
  }

  //_new_h & _inv_new_h
  Util::madd(old_h, d_h, new_h, 3, 3);
  Util::minverse(new_h, inv_new_h, 3);  

  //perturb axes
  for(i = 0; i < 3; i++) 
    Util::vcopy(new_h+3*i, mol.axis(i).x, 3);
  
  //perturb atoms
  for(i = 0; i < natoms; i++) {
    Util::mdotv(inv_old_h, mol.atom(i)->x, s_x, 3, 3);
    Util::mdotv(new_h, s_x, mol.atom(i)->x, 3, 3);
  }
}

void MCAgent::accPertAxes(bool _flag) {}

double MCAgent::urand(double _a, double _b) const {
  return _a+(_b-_a)*Util::RNG::uniform_SFMT();
}

}/* MC */
