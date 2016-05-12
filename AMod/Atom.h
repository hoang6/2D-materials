#ifndef _AMOD_ATOM_H
#define _AMOD_ATOM_H

#include "PTE.h"
#include "../Util/util.h"
#include "../Util/PtrCont.h"
#include <iostream>
#include <vector>
#include <utility>

namespace AMod {

class Atom;
class Arrow;
class Bond;
class BondLength ; 
class BondRange;
typedef Util::PtrCont<Atom> Atoms;
typedef Util::PtrCont<Arrow> Arrows;
typedef Util::PtrCont<Bond> Bonds;
typedef std::pair<int,int> DInt;
typedef std::pair<double,double> DDbl;
typedef std::pair<int,int> IDid;
typedef std::pair<int,bool> IDpos;

class Bond {
public:
  double r; //bond length
  Util::PtrCont<Bond>::ID id;

  Bond();
  Bond(const Bond& bond);
  ~Bond();
  Bond& operator= (const Bond& bond);
  void reset();
  void connect(Atom& atom1, Atom& atom2);
  void connect(Atom& atom1, Atom& atom2, int arid1, int arid2);
  void setup(const double& _r, const double* _dx, const int* _n);
  bool check() const;
  Arrow& arrow(bool pos) const;
  void dumpTxt(std::ostream& fout) const;
  void dumpBin(std::ostream& fout) const;

private:
  Arrow* parrows[2];
};

class Arrow {
public:
  int n[3];
  double dx[3]; //dx = moon()->x-host()->x+n.T
  Util::PtrCont<Arrow>::ID id;

  Arrow();
  Arrow(const Arrow& arrow);
  ~Arrow();
  Arrow& operator= (const Arrow& arrow);
  void reset();
  Arrow& find(int n);
  const Arrow& find(int n) const;
  bool check() const;
  IDid idid() const;
  IDpos idpos() const;
  Bond& bond() const;
  bool position() const;
  Atom& host() const;
  Atom& moon() const;
  Arrow& brother() const;
  void dumpTxt(std::ostream& fout) const;
  void dumpBin(std::ostream& fout) const;

private:
  Atom* patom;
  Bond* pbond; bool pos; //assert(pbond->parrows[pos] == this)
  friend class Bond; 
};

class Atom {
public:
  class Data;
  Arrows arrows;  
  Util::PtrCont<Atom>::ID id;
  
  Atom();
  Atom(Data& data);
  Atom(const Atom& atom);
  ~Atom();
  Atom& operator= (const Atom& atom);
  void reset();
  void attach(Data& data);
  Data& data();
  const Data& data() const;
  Data* operator-> ();
  const Data* operator-> () const;
  int arrowID(int _atomID) const;
  void dumpTxt(std::ostream& fout) const;
  void dumpBin(std::ostream& fout) const;

  class Data {
  public:
    int type;
    double x[3];
    double fpot[3];
    double fext[3];
    int fixed[3];
    Util::PtrCont<Data>::ID id;

    Data();
    Data(int _type, const double* _x);
    Data(int _type, double _x, double _y, double _z);
    Data(const Data& data);
    Data& operator= (const Data& data);
    void reset();
    double ftotal(int i) const;
    double& operator[] (int n);
    const double& operator[] (int n) const;
    void dumpTxt(std::ostream& fout) const;
    void dumpBin(std::ostream& fout) const;
  };

private:
  Data* pdata;
};

class BondLength {
public:
  BondLength();
  BondLength(const BondLength& bondLength);
  ~BondLength();
  BondLength& operator= (const BondLength& bondLength);
  void reset();
  void insert(int _type1, int _type2, double _r);
  void erase(int _type1, int _type2);
  double value(int _type1, int _type2) const;

private:
  double val[PTE::nAtomTypes][PTE::nAtomTypes];
};

class BondRange {//[min, max)
public:
  BondRange();
  BondRange(const BondRange& bondRange);
  ~BondRange();
  BondRange& operator= (const BondRange& bondRange);
  void reset();
  void insert(int _type1, int _type2, double _rmin, double _rmax);
  void erase(int _type1, int _type2);
  double min(int _type1, int _type2) const;
  double max(int _type1, int _type2) const;
  bool ok(int _type1, int _type2, double r) const;

private:
  BondLength lmin;
  BondLength lmax;
};

std::ostream& operator<< (std::ostream& fout, const Bond& bond);

std::ostream& operator<< (std::ostream& fout, const Arrow& arrow);

std::ostream& operator<< (std::ostream& fout, const Atom& atom);

std::ostream& operator << (std::ostream& fout, const Atom::Data& data);

/************************************************************/
inline bool Bond::check() const { return parrows[0]!=NULL && parrows[1]!=NULL; }

inline Arrow& Bond::arrow(bool pos) const { return *parrows[pos]; }

/************************************************************/
inline Arrow& Arrow::find(int n) { return patom->arrows[Util::cycle1D(id()+n, int(patom->arrows.size()))]; }

inline const Arrow& Arrow::find(int n) const { return patom->arrows[Util::cycle1D(id()+n, int(patom->arrows.size()))]; }

inline bool Arrow::check() const { return patom!=NULL && pbond!=NULL && &pbond->arrow(pos)==this; }

inline IDid Arrow::idid() const { return IDid(patom->id(), id()); }

inline IDpos Arrow::idpos() const { return IDpos(pbond->id(), pos); }

inline Bond& Arrow::bond() const { return *pbond; }

inline bool Arrow::position() const { return pos; }

inline Atom& Arrow::host() const { return *patom; }

inline Atom& Arrow::moon() const { return *pbond->arrow(1-pos).patom; }

inline Arrow& Arrow::brother() const { return pbond->arrow(1-pos); }

/************************************************************/
inline void Atom::attach(Data& data) { pdata = &data; }

inline Atom::Data& Atom::data() { return *pdata; }

inline const Atom::Data& Atom::data() const { return *pdata; }

inline Atom::Data* Atom::operator-> () { return pdata; }

inline const Atom::Data* Atom::operator-> () const { return pdata; }


/************************************************************/
inline double Atom::Data::ftotal(int i) const { return fpot[i]+fext[i]; }

inline double& Atom::Data::operator[] (int n) { return x[n]; }

inline const double& Atom::Data::operator[] (int n) const { return x[n]; }

/********** overload << for Bond, Arrow & Atom ***********/
inline std::ostream& operator<< (std::ostream& fout, const Bond& bond) {
  bond.dumpTxt(fout);
  return fout;
}

inline std::ostream& operator<< (std::ostream& fout, const Arrow& arrow) {
  arrow.dumpTxt(fout);
  return fout;
}

inline std::ostream& operator<< (std::ostream& fout, const Atom& atom) {
  atom.dumpTxt(fout);
  return fout;
}
  
inline std::ostream& operator << (std::ostream& fout, const Atom::Data& data) {
  data.dumpTxt(fout);
  return fout;
}

} /* AMod */

#endif
