#ifndef _MD_NVE_H
#define _MD_NVE_H

#include "Data.h"
#include <iostream>
#include <cmath>

namespace MD {

class NVE {
public:
  NVE();
  virtual ~NVE();
  void reset();
  void init(Data& _data);
  void final();
  virtual void step(double _dt);
  //--- the followings are only available after init ---
  double ePot() const { return epot; }
  double eKin() const { return ekin; }
  double eKinMax() const { return ekinMax; }
  int natoms() const { return na; }
  int ndata() const { return nd; }
  Data& data() { return *pdata; }
  AMod::Molecule& molecule() { return *pmol; }
  double& pos(int _i) { return x[_i]; }        //len = ndata()
  double& vel(int _i) { return (*pvel)[_i]; }  //len = ndata()
  double& acc(int _i) { return a[_i]; }        //len = ndata()
  double& mass(int _i) { return m[_i]; }       //len = natoms()
  AMod::KernelBasic& kernel() { return *pkernel; }
  virtual void dumpEnergy(std::ostream& _fout) const;

protected:
  virtual void computeForce();
  void computeEKin();
  void force2acc();

  Data* pdata;
  AMod::Molecule* pmol;
  Data::Array* pvel;
  AMod::KernelBasic* pkernel;
  int na, nd;
  double* x; //Unit: A
  double* a; //Unit: A/fs^2
  double* m; //Unit: eV*fs^2/A^2
  double ekinMax, ekin, epot;
};

}/* MD */

#endif

