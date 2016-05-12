#ifndef _AMOD_MOL_SERVICE_H
#define _AMOD_MOL_SERVICE_H

#include <cstdlib>

namespace AMod {

class Molecule;
template<class MOLECULE> class MolService;
typedef MolService<Molecule> MolServ;
typedef MolService<const Molecule> ConstMolServ;

template<class MOLECULE>
class MolService {
public:
  MolService();
  MolService(MOLECULE& _mol);
  virtual ~MolService();
  virtual void attach(MOLECULE& _mol);
  MOLECULE& molecule() const;

private:
  MOLECULE* pmol;
};

/************************************************************/ 
template<class MOLECULE>
inline MolService<MOLECULE>::MolService() { pmol = (MOLECULE*)NULL; }

template<class MOLECULE>
inline MolService<MOLECULE>::MolService(MOLECULE& _mol) { attach(_mol); }

template<class MOLECULE>
inline MolService<MOLECULE>::~MolService() {}

template<class MOLECULE>
inline void MolService<MOLECULE>::attach(MOLECULE& _mol) { pmol = &_mol; }

template<class MOLECULE>
inline MOLECULE& MolService<MOLECULE>::molecule() const { return *pmol; }

} /* AMod */

#endif
