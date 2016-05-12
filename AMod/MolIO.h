#ifndef _AMOD_MOL_IO_H
#define _AMOD_MOL_IO_H

#include "MolService.h"
#include "../Util/Accessors.h"
#include <iostream>

namespace AMod {

typedef enum {SUCCESS = 0, ISTREAM_ERROR, 
	      BOX_ERROR, ATOM_ERROR, BOND_ERROR} MolIStatus;

class MolI: public MolServ {
public:
  MolI();
  MolI(Molecule& _mol);
  virtual ~MolI();
  void reset();
  static int sourceTxt(std::istream& fin, Molecule& mol);
  static int sourceBin(std::istream& fin, Molecule& mol);
  int sourceTxt(std::istream& fin);
  int sourceBin(std::istream& fin);
};

/************************************************************/
class MolO: public ConstMolServ {
public:
  MolO();
  MolO(const Molecule& _mol);
  virtual ~MolO();
  void reset();
  static void setFormat(std::ostream& fout, int _precision = 9);
  static void dumpTxt(std::ostream& fout, const Molecule&mol, const MolO& molO);
  static void dumpTxt(std::ostream& fout, const Molecule& mol, 
		      bool _dumpCh = true);
  static void dumpBin(std::ostream& fout, const Molecule& mol);
  void dumpTxt(std::ostream& fout) const;
  void dumpBin(std::ostream& fout) const;
  /****** set & get status ******/
  Util::Accessors<bool> dumpChains;
};

/************************************************************/
class MolIO: public MolI, public MolO {
public:
  MolIO();
  MolIO(Molecule& _mol);
  virtual ~MolIO();
  virtual void attach(Molecule& _mol);
  void reset();
  using MolI::molecule;

protected:
  using MolO::attach;
};

/************************************************************/
/********** class MolI**********/
inline MolI::~MolI() {}

inline int MolI::sourceTxt(std::istream& fin) { 
  return sourceTxt(fin, molecule()); }

inline int MolI::sourceBin(std::istream& fin) { 
  return sourceBin(fin, molecule()); }

/********** class MolO**********/
inline MolO::~MolO() {}

inline void MolO::dumpTxt(std::ostream& fout) const { 
  dumpTxt(fout, molecule(), dumpChains()); }

inline void MolO::dumpBin(std::ostream& fout) const { 
  dumpBin(fout, molecule()); }

/********** class MolIO**********/
inline MolIO::~MolIO() {}

} /* AMod */

#endif
