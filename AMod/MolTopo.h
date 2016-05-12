#ifndef _AMOD_MOL_TOPO_H
#define _AMOD_MOL_TOPO_H

#include "MolService.h"
#include "Atom.h"
#include "Trans.h"
#include "../Util/Accessors.h"
#include "../Util/PtrCont.h"
#include "../Util/Array2D.h"
#include <map>
#include <list>
#include <vector>

namespace AMod {

class Molecule;

typedef enum {MODE_NULL = 0, MODE_FREEZE, MODE_CELL} MolTopoMode;

typedef enum {CASE_AXES = 0, CASE_LONELY, CASE_BRUTAL, CASE_GENERAL, 
	      CASE_CHANGE_MODE, CASE_ROTATION, CASE_TRANSLATION } MolTopoCase;

class MolTopo: public MolServ {
public:
  MolTopo();
  MolTopo(Molecule& _mol);
  MolTopo(const MolTopo& molTopo);
  MolTopo& operator= (const MolTopo& molTopo);
  void reset();
  void clear();
  /********** sync & update **********/
  void sync();            //will not change existing topology
  void update(MolTopoCase _case); //update for global changes
  /********** new/del/move atom **********/
  Atom& newAtom(int _atomID, Atom::Data& _atomData);
  void deleteAtom(int _atomID);
  Atom& moveAtom(int _atomID);
  /********** make/new/move lonely atom **********/
  Atom& newLonelyAtom(int _atomID, Atom::Data& _atomData);
  Atom& moveLonelyAtom(int _atomID);
  Atom& makeLonelyAtom(int _atomID);
  /********** new/del bond **********/ 
  Bond& newBond(int _atomID1, int _atomID2);
  Bond& newBond(int _atomID1, int _atomID2, const int* _n);
  void deleteBond(int _bondID);
  /********** metric **********/
  double distance(int _atomID, const double _x[3]) const;
  double distance(int _atomID1, int _atomID2) const;
  void displacement(int _atomID1, int _atomID2, double* _dx) const;
  /********** utilities **********/
  bool bondok(int _bondID) const; //is bondRange().ok?
  int bondID(int _atomID1, int _atomID2) const;
  int atomID(const double* _x) const;
  int atomID(const double* _x, double& _rmin) const;
  void around(const double* _x, int _atomID);
  void around(const double* _x);
  /********** set & get **********/
  int natoms() const;
  int nbonds() const;
  Atom& atom(int _atomID);
  Bond& bond(int _bondID);  
  const Atom& atom(int _atomID) const;
  const Bond& bond(int _bondID) const;
  const Trans& trans() const;

  Util::Accessors<MolTopoMode> mode;
  Util::Accessors<BondRange> bondRange;
  Util::Accessors<double> cellLength;

protected:
  struct Cell;
  struct A2CLink;
  typedef std::list<Atom*> AtomPtrList;
  typedef std::map<Cell,AtomPtrList> C2ALinks;
  typedef Util::PtrCont<A2CLink> A2CLinks;
  typedef Util::Array2D<double,3,3> Axes;
  friend class MolChecker; //for debug

  /********** clone topology **********/
  void clone(const MolTopo& molTopo);
  /********** manage cell **********/
  void setupTrans();
  void setupCells(double _cellLength);
  void clearCells();
  void whatCell(const double* _x, Cell& _cell) const;
  void updateCellsNewAtom(int _atomID);
  void updateCellsDeleteAtom(int _atomID);
  void updateCellsMoveAtom(int _atomID);
  void updateSoftConnection(int _bondID);
  void updateHardConnection(int _atomID, bool _single);
  bool newHardBond(int _atomID1, int _atomID2);
  /********** connect atoms **********/
  void setupAxes();
  void setupConnection();
  void updateConnection(MolTopoMode _mode);
  void updateConnection(MolTopoMode _mode, int _atomID);
  void clearConnection();
  void clearConnection(int _atomID);

  Trans tr, ctr;
  Axes cellAxes;
  int ncell[3];
  C2ALinks c2aLinks;
  A2CLinks a2cLinks;
  Util::Accessors<Atoms> atoms;
  Util::Accessors<Bonds> bonds;

  struct Cell {
    int x[3];
    bool operator== (const Cell& rhs) const;
    bool operator< (const Cell& rhs) const;
  };

  struct A2CLink {
    Cell cell;
    Util::PtrCont<A2CLink>::ID id;
  };
};

/************************************************************/
inline int MolTopo::natoms() const { return atoms().size(); }

inline int MolTopo::nbonds() const { return bonds().size(); }

inline Atom& MolTopo::atom(int _atomID) { return atoms()[_atomID]; }

inline Bond& MolTopo::bond(int _bondID) { return bonds()[_bondID]; }

inline const Atom& MolTopo::atom(int _atomID) const { return atoms()[_atomID]; }

inline const Bond& MolTopo::bond(int _bondID) const { return bonds()[_bondID]; }

inline const Trans& MolTopo::trans() const { return tr; }

}/* AMod */

#endif
