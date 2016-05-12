#ifndef _MD_NVE_1HE_H
#define _MD_NVE_1HE_H

#include "NVE.h"
#include "LJ.h"

namespace MD {

class NVE1He: public NVE {
public:
  struct Para;
  NVE1He();
  virtual ~NVE1He();
  void reset();
  //Helium will be appended as the last atom of mol
  void init(Data& _data, const Para& _para);
  void final();
  double eHeC() const { return ehec; }
  AMod::Atom& atomHe() { return molecule().atom(atomHeID); };
  double& posHe(int _i) { return NVE::pos(ndata()-3+_i); }
  double& velHe(int _i) { return NVE::vel(ndata()-3+_i); }
  double& accHe(int _i) { return NVE::acc(ndata()-3+_i); }
  double& massHe() { return mass(natoms()-1); }
  virtual void dumpEnergy(std::ostream& _fout) const;

protected:
  virtual void computeForce();

  const LJCutoff ljc;
  AMod::MolTopo topoHeC;
  double ehec;
  int atomHeID;
};

struct NVE1He::Para {
  double posHe[3]; double velHe[3];
  //------------------------------------------------------------
  Para() { reset(); }
  Para(const double* _posHe, const double* _velHe) {
    for(int i = 0; i < 3; i++) {posHe[i] = _posHe[i]; velHe[i] = _velHe[i];}}
  void reset() {for(int i = 0; i < 3; i++) {posHe[i] = 0.0; velHe[i] = 0.0;}}
};

}/* MD */

#endif
