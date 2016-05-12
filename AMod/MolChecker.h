#ifndef _AMOD_MOL_CHECKER
#define _AMOD_MOL_CHECKER

#include "MolService.h"

namespace AMod {

class MolChecker: public ConstMolServ {
public:
  MolChecker();
  MolChecker(const Molecule& _mol);
  void reset();
  bool neighbCheck() const;
  bool connectionCheck() const;
  bool topoCheck() const;
  bool check() const;
  void printAtom(int _atomID) const;

protected:
  template <class T>
  static void printVec(const T* v, long len);
  void showMessage(const char* whatCheck, bool passed) const;
};

} /* AMod */

#endif
