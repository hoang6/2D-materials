#include "MolTopo.h"
#include "Molecule.h"
#include "../Util/util.h"
#include <iostream>
#include <set>
#include <cmath>

namespace AMod {

MolTopo::MolTopo() { reset(); }

MolTopo::MolTopo(Molecule& _mol): MolServ(_mol) { reset(); }

MolTopo::MolTopo(const MolTopo& molTopo) { *this = molTopo; }

MolTopo& MolTopo::operator= (const MolTopo& molTopo) {
  if(this != &molTopo) {
    MolServ::operator=(molTopo);
    mode = molTopo.mode;
    bondRange = molTopo.bondRange;
    cellLength = molTopo.cellLength;
    clone(molTopo);
  }
  return *this;
}

void MolTopo::reset() {
  mode(MODE_CELL);
  bondRange().reset();
  bondRange().insert(PTE::Carbon,   PTE::Carbon,   0.0, 2.0);
  bondRange().insert(PTE::Carbon,   PTE::Hydrogen, 0.0, 1.8);
  bondRange().insert(PTE::Hydrogen, PTE::Hydrogen, 0.0, 1.7);
  cellLength(2.05);
  clear();
}

void MolTopo::clear() {
  tr.reset();
  ctr.reset();
  cellAxes = Axes(1.0,0.0);
  ncell[0] = ncell[1] = ncell[2] = -1;
  clearCells();
  atoms().clear();
  bonds().clear();
}

void MolTopo::sync() {
  Molecule& mol = molecule();
  int k, kmax = mol.atomsData().size();
  atoms().resize(kmax);
  for(k = 0; k < kmax; k++) atom(k).attach(mol.atomsData()[k]);
}

void MolTopo::update(MolTopoCase _case) {
  if(_case == CASE_AXES) {
    setupAxes();
  }
  else if(_case == CASE_LONELY) {
    clearConnection();
  }
  else if(_case == CASE_BRUTAL) {
    if(mode() == MODE_FREEZE)
      setupConnection();
    else
      updateConnection(mode());
  }
  else if(_case == CASE_GENERAL) {
    updateConnection(mode());
  }
  else if(_case == CASE_CHANGE_MODE) {
    if(mode() == MODE_CELL)
      setupCells(cellLength());
    else 
      clearCells();
  }
  else if(_case == CASE_ROTATION) {
    setupAxes();
    updateConnection(MODE_FREEZE);
  }
  else if(_case == CASE_TRANSLATION) {
    if(mode() == MODE_CELL)
      setupCells(cellLength());
  }
}

Atom& MolTopo::newAtom(int _atomID, Atom::Data& _atomData) {
  Atom& a = (_atomID<0 ? atoms().insert() : atoms().insert(_atomID));
  a.reset();
  a.attach(_atomData);
  if(mode() == MODE_NULL) {
    for(int k = 0, kmax = natoms(); k < kmax; k++) {
      if(a.id() != k) newHardBond(a.id(), k);
    }
  }
  else if(mode() == MODE_CELL) {
    updateCellsNewAtom(a.id());
    updateHardConnection(a.id(), true);
  }
  return a;
}

void MolTopo::deleteAtom(int _atomID) {
  if(mode() == MODE_CELL)
    updateCellsDeleteAtom(_atomID);
  clearConnection(_atomID);
  atoms().erase(_atomID);
}

Atom& MolTopo::moveAtom(int _atomID) {
  updateConnection(mode(), _atomID);
  return atom(_atomID);
}

Atom& MolTopo::newLonelyAtom(int _atomID, Atom::Data& _atomData) {
  Atom& a = (_atomID<0 ? atoms().insert() : atoms().insert(_atomID));
  a.reset();
  a.attach(_atomData);
  if(mode() == MODE_CELL) updateCellsNewAtom(a.id());
  return a;
}

Atom& MolTopo::moveLonelyAtom(int _atomID) {
  if(mode() == MODE_CELL)
    updateCellsMoveAtom(_atomID);
  return atom(_atomID);
}

Atom& MolTopo::makeLonelyAtom(int _atomID) {
  clearConnection(_atomID);
  return atom(_atomID);
}

Bond& MolTopo::newBond(int _atomID1, int _atomID2) {
  Bond& b = bonds().insert();
  b.connect(atom(_atomID1), atom(_atomID2));
  updateSoftConnection(b.id());
  return b;
}

Bond& MolTopo::newBond(int _atomID1, int _atomID2, const int* _n) {
  Molecule& mol = molecule();
  Atom& a1 = atom(_atomID1);
  Atom& a2 = atom(_atomID2);
  Bond& b = bonds().insert();
  double r, dx[3];
  b.connect(a1, a2);
  Util::vsub(a2->x, a1->x, dx, 3);
  for(int k = 0; k < mol.naxes(); k++) {
    Util::kdotvadd(1.0, dx, double(_n[k]), mol.axis(k).x, dx, 3);
  }
  r = Util::norm2(dx, 3);
  b.setup(r, dx, _n);
  return b;
}

void MolTopo::deleteBond(int _bondID) {
  Bond& b = bond(_bondID);
  for(int k = 0; k < 2; k++)
    b.arrow(k).host().arrows.erase(b.arrow(k).id());
  bonds().erase(_bondID);
}

double MolTopo::distance(int _atomID, const double _x[3]) const {
  const double* xa = atom(_atomID)->x;
  double xc[3] = {_x[0], _x[1], _x[2]};
  tr.findImage(xc, xa);
  return std::sqrt(tr.distSquare());
}

double MolTopo::distance(int _atomID1, int _atomID2) const {
  const double* x1 = atom(_atomID1)->x;
  const double* x2 = atom(_atomID2)->x;
  double cx2[3] = {x2[0], x2[1], x2[2]};
  tr.findImage(cx2, x1);
  return std::sqrt(tr.distSquare());
}

void MolTopo::displacement(int _atomID1, int _atomID2, double* _dx) const {
  const double* x1 = atom(_atomID1)->x;
  const double* x2 = atom(_atomID2)->x;
  double cx2[3] = {x2[0], x2[1], x2[2]};
  tr.findImage(cx2, x1);
  for(int k = 0; k < 3; k++) _dx[k] = tr.displacement(k);
}

bool MolTopo::bondok(int _bondID) const {
  const Bond& b = bond(_bondID);
  return bondRange().ok(b.arrow(0).host()->type,b.arrow(1).host()->type,b.r);
}

int MolTopo::bondID(int _atomID1, int _atomID2) const {
  const Atom& a1 = atom(_atomID1);
  int _arrowID = a1.arrowID(_atomID2);
  int _bondID = -1;
  if(_arrowID >= 0) _bondID = a1.arrows[_arrowID].bond().id();
  return _bondID;
}

int MolTopo::atomID(const double* _x) const {
  double _rmin;
  return atomID(_x, _rmin);
}

int MolTopo::atomID(const double* _x, double& _rmin) const {
  int _atomID = -1;
  double _r2, _r2min = 0.0;
  for(int k = 0, kmax = natoms(); k < kmax; k++) {
    double _cx[3] = {_x[0], _x[1], _x[2]};
    tr.findImage(_cx, atom(k)->x);
    _r2 = tr.distSquare();
    if(k==0 || _r2 < _r2min) {
      _atomID = k;
      _r2min = _r2;
    }
  }
  _rmin = std::sqrt(_r2min);
  return _atomID;
}

void MolTopo::around(const double* _x, int _atomID) {
  tr.findImage(atom(_atomID)->x, _x);
  moveAtom(_atomID);
}

void MolTopo::around(const double* _x) {
  for(int k = 0, kmax = natoms(); k < kmax; k++) around(_x, k);
}

void MolTopo::clone(const MolTopo& molTopo) {
  /* Atoms */
  atoms = molTopo.atoms;
  /* Bonds */
  int k, kmax;
  IDid idid[2];
  bonds = molTopo.bonds;
  kmax = bonds().size();
  for(k = 0; k < kmax; k++) {
    idid[0] = molTopo.bond(k).arrow(0).idid();
    idid[1] = molTopo.bond(k).arrow(1).idid();
    bond(k).connect(atom(idid[0].first), atom(idid[1].first),
		    idid[0].second, idid[1].second); 
  }
  /* Cell */
  C2ALinks::iterator it;
  AtomPtrList::iterator ait;
  tr = molTopo.tr;
  ctr = molTopo.ctr;
  cellAxes = molTopo.cellAxes;
  for(k = 0; k < 3; k++) 
    ncell[k] = molTopo.ncell[k];
  c2aLinks = molTopo.c2aLinks;
  a2cLinks = molTopo.a2cLinks;
  for(it = c2aLinks.begin(); it != c2aLinks.end(); it++) {
    AtomPtrList& apl = (*it).second;
    for(ait = apl.begin(); ait != apl.end(); ait++) {
      (*ait) = &atom((*ait)->id());
    }
  }
}

void MolTopo::setupTrans() {
  tr.init(molecule().naxes(), molecule().axes()); 
}

void MolTopo::setupCells(double _cellLength) {
  const Molecule& mol = molecule();
  int j, k, kmax;
  double axesLengths[3];
  Axes dualCellAxes;
  double temp;
  Cell ctemp;
  //setup cellAxes
  for(k = 0, kmax = mol.naxes(); k < kmax; k++)
    for(j = 0; j < 3; j++)
      cellAxes[k][j] = mol.axis(k).x[j];
  Util::build_axes(kmax, cellAxes);
  //find axesLengths of cellAxes
  for(k = 0; k < 3; k++) {
    Util::normalize(cellAxes[k], 3);
  }
  Util::dual_axes(cellAxes[0], cellAxes[1], cellAxes[2],
		  dualCellAxes[0], dualCellAxes[1], dualCellAxes[2]);
  for(k = 0; k < 3; k++) {
    ncell[k] = -1; //no periodicity in this direction
    axesLengths[k] = _cellLength*Util::norm2(dualCellAxes[k], 3);
  }
  for(k = 0, kmax = mol.naxes(); k < kmax; k++) {
    temp = Util::norm2(mol.axis(k).x, 3);
    ncell[k] = int(temp/axesLengths[k]);
    axesLengths[k] = temp/ncell[k];
  }
  //scale cellAxes
  for(k = 0; k < 3; k++)
    Util::kdotv(axesLengths[k], cellAxes[k], cellAxes[k], 3);
  //fill c2aLinks and a2cLinks
  clearCells();
  a2cLinks.resize(natoms());
  ctr.init(3, cellAxes); //for whatCell(double*, Cell&)
  for(k = 0, kmax = natoms(); k < kmax; k++) {
    whatCell(atom(k)->x,ctemp);
    c2aLinks[ctemp].push_back(&atom(k));
    a2cLinks[k].cell = ctemp;
  }
}

void MolTopo::clearCells() {
  c2aLinks.clear();
  a2cLinks.clear();
}

void MolTopo::whatCell(const double* _x, Cell& _cell) const {
  int i;
  double tx[3];
  ctr.proj(_x, tx);
  for(i = 0; i < 3; i++)
    _cell.x[i] = int(std::floor(tx[i]));
  for(i = 0; i < 3; i++)
    if(ncell[i] > 0) {
      _cell.x[i] = Util::cycle1D(_cell.x[i], ncell[i]);
    }
}

void MolTopo::updateCellsNewAtom(int _atomID) {
  Atom& a = atom(_atomID);
  Cell& c = a2cLinks.insert(_atomID).cell;
  whatCell(a->x, c);
  c2aLinks[c].push_back(&a);
}

void MolTopo::updateCellsDeleteAtom(int _atomID) {
  Atom& a = atom(_atomID);
  Cell c = a2cLinks[_atomID].cell;
  C2ALinks::iterator it = c2aLinks.find(c);
  (*it).second.remove(&a);
  if((*it).second.empty()) c2aLinks.erase(it);
  a2cLinks.erase(_atomID);
}

void MolTopo::updateCellsMoveAtom(int _atomID) {
  Atom& a = atom(_atomID);
  Cell oldc, newc;
  oldc = a2cLinks[_atomID].cell;
  whatCell(a->x, newc);
  if(oldc == newc) return;
  C2ALinks::iterator it = c2aLinks.find(oldc);
  (*it).second.remove(&a);
  if((*it).second.empty()) c2aLinks.erase(it);
  c2aLinks[newc].push_back(&a);
  a2cLinks[_atomID].cell = newc;
}

void MolTopo::updateSoftConnection(int _bondID) {
  int k;
  double cx2[3];
  const double *x1, *x2;
  Bond& b = bond(_bondID);
  x1 = b.arrow(0).host()->x;
  x2 = b.arrow(1).host()->x;
  for(k = 0; k < 3; k++) cx2[k] = x2[k];
  tr.findImage(cx2, x1);
  //set Arrow::dx & Arrow:dn
  b.setup(std::sqrt(tr.distSquare()), tr.displacement(), tr.nT());
}

void MolTopo::updateHardConnection(int _atomID, bool _single) {
  int i, j, k;
  int _atomID2;
  C2ALinks::iterator cit;
  AtomPtrList::iterator ait;
  Cell c, ct = a2cLinks[_atomID].cell;
  std::set<Cell> cset;
  std::set<Cell>::iterator sit;
  //new bonds (assume atom is in the right cell)
  for(i = -1; i <= +1; i++)
    for(j = -1; j <= +1; j++)
      for(k = -1; k <= +1; k++) {
	c = ct;
	c.x[0] += i;
	c.x[1] += j;
	c.x[2] += k;
	if(ncell[0] > 0) c.x[0] = (c.x[0]+ncell[0])%ncell[0];
	if(ncell[1] > 0) c.x[1] = (c.x[1]+ncell[1])%ncell[1];
	if(ncell[2] > 0) c.x[2] = (c.x[2]+ncell[2])%ncell[2];
	cset.insert(c);
      }
  for(sit = cset.begin(); sit != cset.end(); sit++) {
    cit = c2aLinks.find(*sit);
    if(cit != c2aLinks.end()) {
      AtomPtrList& apl = (*cit).second;
      for(ait = apl.begin(); ait != apl.end(); ait++) {
	_atomID2 = (*ait)->id();
	if((_single&&(_atomID!=_atomID2)) || (_atomID<_atomID2))
	  newHardBond(_atomID, _atomID2);
      }
    }
  }
}

bool MolTopo::newHardBond(int _atomID1, int _atomID2) {
  int k;
  double *x1, x2[3], rmin;
  Atom& a1 = atom(_atomID1);
  Atom& a2 = atom(_atomID2);
  x1 = a1->x;
  for(k = 0; k < 3; k++) x2[k] = a2->x[k];
  if(tr.findImage(x2, x1, bondRange().max(a2->type,a1->type))) {
    rmin = bondRange().min(a2->type,a1->type);
    if(tr.distSquare() >= rmin*rmin) {
      //set the geometry part of a bond
      Bond& bond = bonds().insert();
      bond.connect(a1, a2);
      //set Arrow::dx & Arrow:dn
      bond.setup(std::sqrt(tr.distSquare()), tr.displacement(), tr.nT());
      return true;
    }
  }//end if(trans.findImage...
  return false;
}

void MolTopo::setupAxes() {
  setupTrans();
  if(mode() == MODE_CELL) 
    setupCells(cellLength());
}

void MolTopo::setupConnection() {
  int k1, k2, kmax = natoms();
  clearConnection();
  for(k1 = 0; k1 < kmax; k1++)
    for(k2 = k1+1; k2 < kmax; k2++) 
      newHardBond(k1, k2);
}

void MolTopo::updateConnection(MolTopoMode _mode) {
  int k, kmax;
  if(_mode == MODE_NULL) {
    setupConnection();
  }
  else if(_mode == MODE_CELL) {
    clearConnection();
    kmax = natoms();
    for(k = 0; k < kmax; k++) updateCellsMoveAtom(k);
    for(k = 0; k < kmax; k++) updateHardConnection(k, false);
  }
  else if(_mode == MODE_FREEZE) {
    kmax = nbonds();
    for(k = 0; k < kmax; k++) updateSoftConnection(k);
  }
}

void MolTopo::updateConnection(MolTopoMode _mode, int _atomID) {
  Atom& a = atom(_atomID);
  if(_mode == MODE_NULL) {
    clearConnection(_atomID);
    for(int k = 0, kmax = natoms(); k < kmax; k++) {
      if(a.id() != k) newHardBond(a.id(), k);
    }
  }
  else if(_mode == MODE_CELL){
    clearConnection(_atomID);
    updateCellsMoveAtom(_atomID);
    updateHardConnection(_atomID, true);
  }
  else if(_mode == MODE_FREEZE) {
    for(int k = 0, kmax = a.arrows.size(); k < kmax; k++)
      updateSoftConnection(a.arrows[k].bond().id());
  }
}

void MolTopo::clearConnection() {
  for(int k = 0, kmax = natoms(); k < kmax; k++) 
    atom(k).arrows.clear();
  bonds().clear();  
}

void MolTopo::clearConnection(int _atomID) {
  Atom& a = atom(_atomID);
  int i, imax = a.arrows.size();
  if(imax == 0) return;
  for(i = 0; i < imax; i++)
    deleteBond(a.arrows[0].bond().id());
  a.arrows.clear();
}

bool MolTopo::Cell::operator== (const Cell& rhs) const {
  for(int k = 0; k < 3; k++)
    if(x[k] != rhs.x[k]) return false;
  return true;
}

bool MolTopo::Cell::operator< (const MolTopo::Cell& rhs) const {
  for(int k = 0; k < 3; k++) {
    if(x[k] < rhs.x[k]) return true;
    if(x[k] > rhs.x[k]) return false;
  }
  return false;
}

}/* AMod */
