#include "MolScaler.h"
#include "Molecule.h"
#include "../Util/util.h"

namespace AMod {

MolScaler::MolScaler() { reset(); }

MolScaler::MolScaler(Molecule& _mol): ConstMolServ(_mol) { init(); }
 
void MolScaler::reset() {
  Util::vassign(0.0, h, 9);
  Util::vassign(0.0, inv_h, 9);
}

void MolScaler::attach(const Molecule& _mol) {
  MolService<const Molecule>::attach(_mol);
  init();
}

void MolScaler::init() {
  for(int k = 0; k < 3; k++) 
    Util::vcopy(molecule().axis(k).x, h+3*k, 3);
  Util::minverse(h, inv_h, 3);
}

void MolScaler::scale(int _atomID, double _sx[3]) const {
  Util::mdotv(inv_h, molecule().atom(_atomID)->x, _sx, 3, 3);
}

}/* AMod */
