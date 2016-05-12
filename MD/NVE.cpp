#include "NVE.h"
#include "../AMod/KernelREBO.h"
#include "../AMod/KernelLCBOPIIN.h"

namespace MD {

NVE::NVE() { reset(); }

NVE::~NVE() {}

void NVE::reset() {
  pdata = NULL;
  pmol = NULL;
  pvel = NULL;
  pkernel = NULL;
  na = 0; nd = 0;
  x = NULL; a = NULL; m = NULL;
  ekin = 0.0; epot = 0.0;
}

void NVE::init(Data& _data) {
  _data.forceConsistency();
  pdata = &_data;
  pmol = &_data.mol;
  pvel = &_data.vel;
  if(_data.pot == AMod::REBO)
    pkernel = new AMod::KernelREBO(*pmol);
  else if(_data.pot == AMod::LCBOPIIN)
    pkernel = new AMod::KernelLCBOPIIN(*pmol);
  pkernel->init();

  na = pmol->natoms();
  nd = 3*na;

  x = new double[nd];
  for(int i = 0; i < nd; i++)
    x[i] = pmol->atom(i/3)->x[i%3];

  m = new double[na];
  for(int i = 0; i < na; i++) m[i] = _data.mass(i);

  a = new double[nd];

  computeEKin();
  computeForce();
  force2acc();
}

void NVE::final() {
  delete [] a;
  delete [] m;
  delete [] x;
  pkernel->final();
  delete pkernel;
}

void NVE::step(double _dt) {
  for(int i = 0; i < ndata(); i++) { 
    vel(i) += acc(i)*(_dt/2.0);
    pos(i) += vel(i)*_dt;
  }
  computeForce();
  force2acc();
  for(int i = 0; i < ndata(); i++)
    vel(i) += acc(i)*(_dt/2.0);
  computeEKin();
}

void NVE::dumpEnergy(std::ostream& _fout) const {
  _fout << std::showpos << eKin() << "\t" << ePot() << std::endl;
}

void NVE::computeForce() {
  pkernel->update(x);
  pkernel->energy();
  pkernel->denergy();
  epot = molecule().potential();
}

void NVE::computeEKin() {
  int i;
  double tmp;
  ekin = 0.0, ekinMax = 0.0;
  for(i = 0; i < ndata(); i++) {
    tmp = 0.5*mass(i/3)*vel(i)*vel(i);
    ekin += tmp; if(tmp > ekinMax) ekinMax = tmp;
  }
}

void NVE::force2acc() {
  for(int i = 0; i < natoms(); i++)
    for(int j = 0; j < 3; j++)
      acc(3*i+j) = molecule().atom(i)->fpot[j]/m[i];
}

}/* MD */
