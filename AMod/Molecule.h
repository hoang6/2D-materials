#ifndef _AMOD_MOLECULE_H
#define _AMOD_MOLECULE_H

#include "Atom.h"
#include "MolTopo.h"
#include "MolIO.h"
#include "../Util/Accessors.h"
#include <vector>
#include <string>
#include <iostream>

namespace AMod {

class Molecule;
typedef Util::PtrCont<Molecule> Molecules;

class Molecule {
public:
  typedef Atom::Data Axis;
  typedef Atom::Data AData;
  typedef Util::PtrCont<Axis> Axes;
  typedef Util::PtrCont<AData> AtomsData;

  Molecule();
  Molecule(const Molecule& mol);
  Molecule& operator= (const Molecule& mol);
  void reset();
  void clear();
  /********** new/del/move atoms **********/
  Atom& newAtom(int _atomID, const AData& _atomData);
  void deleteAtom(int _atomID);
  Atom& moveAtom(int _atomID, const double* _x);
  /********** make/new/move lonely atoms **********/
  Atom& newLonelyAtom(int _atomID, const AData& _atomData);
  Atom& moveLonelyAtom(int _atomID, const double* _x);
  Atom& makeLonelyAtom(int _atomID);
  /********** utilities **********/
  void setStdAxes();
  void setFixedAxes(int _dim);
  int countAtoms(int _atomType) const;
  /********** set & get **********/
  int naxes() const;
  int natoms() const;
  int nbonds() const;
  Axis& axis(int _axisID);
  Atom& atom(int _atomID);
  Bond& bond(int _bondID);
  const Axis& axis(int _axisID) const;
  const Atom& atom(int _atomID) const;
  const Bond& bond(int _bondID) const;
  double totalEnergy() const;
  void axes2array(double* _h) const;
  double volume() const; //assume naxes() == 3

  Util::Accessors<Axes> axes;
  Util::Accessors<AtomsData> atomsData;
  Util::Accessors<double> potential;
  Util::Accessors<double> work;
  Util::Accessors<MolIO> io;
  Util::Accessors<MolTopo> topo;
  Util::PtrCont<Molecule>::ID id;

private:
  void attach();
  void getRange(const double _n[], double& _rmin, double& _rmax) const;
};

std::ostream& operator<< (std::ostream& _fout, const Molecule& _mol);

std::istream& operator>> (std::istream& _fin, Molecule& _mol);

std::ostream& operator<< (std::ostream& _fout, const Molecules& _mols);

std::istream& operator>> (std::istream& _fin, Molecules& _mols);

/************************************************************/
inline int Molecule::naxes() const { return axes().size(); }

inline int Molecule::natoms() const { return atomsData().size(); }
  
inline int Molecule::nbonds() const { return topo().nbonds(); }
  
inline Molecule::Axis& Molecule::axis(int _axisID) { return axes()[_axisID]; }
 
inline Atom& Molecule::atom(int _atomID) { return topo().atom(_atomID); }

inline Bond& Molecule::bond(int _bondID) { return topo().bond(_bondID); }

inline const Molecule::Axis& Molecule::axis(int _axisID) const { 
  return axes()[_axisID]; }

inline const Atom& Molecule::atom(int _atomID) const { 
  return topo().atom(_atomID); }

inline const Bond& Molecule::bond(int _bondID) const { 
  return topo().bond(_bondID); }

inline double Molecule::totalEnergy() const { return potential()+work(); }

inline std::ostream& operator<< (std::ostream& _fout, const Molecule& _mol) {
  _mol.io().dumpTxt(_fout); return _fout;
}

inline std::istream& operator>> (std::istream& _fin, Molecule& _mol) {
  _mol.io().sourceTxt(_fin); return _fin;
}

inline std::ostream& operator<< (std::ostream& _fout, const Molecules& _mols) {
  int k, kmax = _mols.size();
  for(k = 0; k < kmax; k++) _fout << _mols[k]; return _fout;
}

inline std::istream& operator>> (std::istream& _fin, Molecules& _mols) {
  _mols.push_back();
  while(!_mols.back().io().sourceTxt(_fin)) _mols.push_back(); 
  _mols.pop_back(); return _fin;
}

}/* AMod */

#endif
