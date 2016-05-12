#include "Atom.h"
#include "../Util/constants.h"
#include <cmath>

namespace AMod {

/********** class Bond ***********/
Bond::Bond() { 
  reset(); 
  parrows[0] = parrows[1] = NULL;
}

Bond::Bond(const Bond& bond) { *this = bond; }

Bond::~Bond() {}

Bond& Bond::operator= (const Bond& bond) {
  if(this != &bond) {
    r = bond.r;
  }
  return *this;
}

void Bond::reset() { r = 0.0; }

void Bond::connect(Atom& atom1, Atom& atom2) {
  Arrow& ar1 = atom1.arrows.insert();
  Arrow& ar2 = atom2.arrows.insert();
  //bond
  parrows[0] = &ar1;
  parrows[1] = &ar2;
  //ar1
  ar1.patom = &atom1;
  ar1.pbond = this;
  ar1.pos = 0;
  //ar2
  ar2.patom = &atom2;
  ar2.pbond = this;
  ar2.pos = 1;
}

void Bond::connect(Atom& atom1, Atom& atom2, int arid1, int arid2) {
  Arrow& ar1 = atom1.arrows[arid1];
  Arrow& ar2 = atom2.arrows[arid2];
  //bond
  parrows[0] = &ar1;
  parrows[1] = &ar2;
  //ar1
  ar1.patom = &atom1;
  ar1.pbond = this;
  ar1.pos = 0;
  //ar2
  ar2.patom = &atom2;
  ar2.pbond = this;
  ar2.pos = 1;
}

void Bond::setup(const double& _r, const double* _dx, const int* _n) {
  r = _r;
  for(int k = 0; k < 3; k++) {
    parrows[0]->dx[k] = _dx[k];
    parrows[0]->n[k] = _n[k];
    parrows[1]->dx[k] = -_dx[k];
    parrows[1]->n[k] = -_n[k];
  }
}

void Bond::dumpTxt(std::ostream& fout) const {
  if(!check()) return;
  int* n = parrows[0]->n;
  fout << id() << " "
       << parrows[0]->host().id() << " "
       << parrows[1]->host().id() << " "
       << n[0] << " " << n[1] << " " << n[2];
}

void Bond::dumpBin(std::ostream& fout) const {
  if(!check()) return;
  int data[3] = {id(), parrows[0]->host().id(), parrows[1]->host().id()};
  fout.write((char*)data, sizeof(data));
  fout.write((char*)parrows[0]->n, sizeof(parrows[0]->n));
}

/********** class Arrow **********/
Arrow::Arrow() { 
  reset(); 
  patom = NULL;
  pbond = NULL;
  pos = 0;
}

Arrow::Arrow(const Arrow& arrow) { *this = arrow; }

Arrow::~Arrow() {}

Arrow& Arrow::operator= (const Arrow& arrow) {
  if(this != &arrow) {
    for(int j = 0; j < 3; j++) {
      n[j] = arrow.n[j];
      dx[j] = arrow.dx[j];
    }
  }
  return *this;
}

void Arrow::reset() {
  for(int j = 0; j < 3; j++) {
    n[j] = 0;
    dx[j] = 0.0;
  }
}

void Arrow::dumpTxt(std::ostream& fout) const {
  if(!check()) return;
  fout << id() << " "
       << patom->id() << " "
       << pbond->id() << " "
       << pos << " "
       << n[0] << " " << n[1] << " " << n[2] << " "
       << dx[0] << " " << dx[1] << " " << dx[2];
}

void Arrow::dumpBin(std::ostream& fout) const {
  if(!check()) return;
  int data[3] = {id(), patom->id(), pbond->id()};
  fout.write((char*)data, sizeof(data));
  fout.write((char*)pos, sizeof(pos));
  fout.write((char*)n, sizeof(n));
  fout.write((char*)dx, sizeof(dx));
}

/********** class Atom ***********/
Atom::Atom() { reset(); }

Atom::Atom(Data& data) { reset(); attach(data); }

Atom::Atom(const Atom& atom) { *this = atom; }

Atom::~Atom() {}

Atom& Atom::operator= (const Atom& atom) {
  if(this != &atom) {
    arrows = atom.arrows;
  }
  return *this;
}

void Atom::reset() { 
  arrows.clear();
  pdata = NULL;
}

int Atom::arrowID(int _atomID) const {
  int _arrowID = -1;
  for(int k = 0, kmax = arrows.size(); k < kmax; k++) {
    if(arrows[k].moon().id() == _atomID) {
      _arrowID = k;
      break;
    }
  }
  return _arrowID;
}

void Atom::dumpTxt(std::ostream& fout) const {
  data().dumpTxt(fout);
}

void Atom::dumpBin(std::ostream& fout) const {
  data().dumpBin(fout);
}

/********** class Atom::Data ***********/
Atom::Data::Data() { reset(); }

Atom::Data::Data(int _type, const double* _x) {
  reset();
  type = _type;
  for(int k = 0; k < 3; k++) x[k] = _x[k];
}

Atom::Data::Data(int _type, double _x, double _y, double _z) {
  reset();
  type = _type;
  x[0] = _x;
  x[1] = _y;
  x[2] = _z;
}

Atom::Data::Data(const Data& data) { *this = data; }

Atom::Data& Atom::Data::operator= (const Data& data) {
  if(this != &data) {
    type = data.type;
    for(int k = 0; k < 3; k++) {
      x[k] = data.x[k];
      fpot[k] = data.fpot[k];
      fext[k] = data.fext[k];
      fixed[k] = data.fixed[k];
    }
  }
  return *this;
}

void Atom::Data::reset() {
  type = PTE::Unknown;
  for(int k = 0; k < 3; k++) {
    x[k] = 0.0;
    fpot[k] = 0.0;
    fext[k] = 0.0;
    fixed[k] = 0;
  }
}

void Atom::Data::dumpTxt(std::ostream& fout) const {
  fout << id() << " " << type << " "
       << x[0] << " " << x[1] << " " << x[2] << " "
       << fpot[0] << " " << fpot[1] << " " << fpot[2] << " "
       << fext[0] << " " << fext[1] << " " << fext[2] << " "
       << fixed[0] << " " << fixed[1] << " " << fixed[2];
}

void Atom::Data::dumpBin(std::ostream& fout) const {
  int myid = id();
  fout.write((char*)&myid, sizeof(myid));
  fout.write((char*)&type, sizeof(type));
  fout.write((char*)x, sizeof(x));
  fout.write((char*)fpot, sizeof(fpot));
  fout.write((char*)fext, sizeof(fext));
  fout.write((char*)fixed, sizeof(fixed));
}

/********** class BondLength **********/
BondLength::BondLength() { reset(); }

BondLength::BondLength(const BondLength& bondLength) { *this = bondLength; }

BondLength::~BondLength() {}

BondLength& BondLength::operator= (const BondLength& bondLength) {
  if(this != &bondLength) {
    int i, j, n = PTE::nAtomTypes;
    for(i = 0; i < n; i++)
      for(j = 0; j < n; j++) val[i][j] = bondLength.val[i][j];
  }
  return *this;
}

void BondLength::reset() {
  int i, j, n = PTE::nAtomTypes;
  for(i = 0; i < n; i++)
    for(j = 0; j < n; j++) val[i][j] = -1.0;
}

void BondLength::insert(int _type1, int _type2, double _r) {
  val[_type1][_type2] = _r;
  if(_type1 != _type2)
    val[_type2][_type1] = _r;
}

void BondLength::erase(int _type1, int _type2) {
  val[_type1][_type2] = -1.0;
  if(_type1 != _type2)
    val[_type2][_type1] = -1.0;
}

double BondLength::value (int _type1, int _type2) const {
  return val[_type1][_type2];
}

/********** class BondRules **********/
BondRange::BondRange() { reset(); }

BondRange::BondRange(const BondRange& bondRange) { 
  *this = bondRange; 
}

BondRange::~BondRange() {}

BondRange& BondRange::operator= (const BondRange& bondRange) {
  if(this != &bondRange) {
    lmin = bondRange.lmin;
    lmax = bondRange.lmax;
  }
  return *this;
}

void BondRange::reset() {
  lmin.reset();
  lmax.reset();
}

void BondRange::insert(int _type1, int _type2, double _rmin, double _rmax) {
  lmin.insert(_type1, _type2, _rmin);
  lmax.insert(_type1, _type2, _rmax);
}

void BondRange::erase(int _type1, int _type2) {
  lmin.erase(_type1, _type2);
  lmax.erase(_type1, _type2);
}

double BondRange::min(int _type1, int _type2) const {
  double d = lmin.value(_type1, _type2);
  return (d<0?Util::MAX_DOUBLE:d);
}

double BondRange::max(int _type1, int _type2) const {
  double d = lmax.value(_type1, _type2);
  return (d<0?0.0:d);
}

bool BondRange::ok(int _type1, int _type2, double r) const {
  if((r >= min(_type1,_type2)) &&
     (r <  max(_type1,_type2))) return true;
  return false;
}

} /* AMod */
