#include "../AMod/AMod.h"
#include "../MC/MC_LCBOPII.h"
#include "../Util/ODErkck.h"
#include "../Util/util.h"
#include <string>
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
using namespace std;

class HeCScatt {
public:
  typedef Util::ODErkck<12> ODErkck;

  void orbit(const AMod::Molecule& _mol, int _atomID, 
	     double _Ek, double _s, double _q, ODErkck::Soln& _soln) {
    init(_mol, _atomID);
    const double* x = mol.atom(atomID)->x;

    ODErkck::Data u0;
    u0[0] = x[0]+_s*cos(_q);
    u0[1] = x[1]+_s*sin(_q);
    u0[2] = x[2]+dFar;
    u0[3] = x[0];
    u0[4] = x[1];
    u0[5] = x[2];
    u0[6] = 0.0;
    u0[7] = 0.0;
    u0[8] = -sqrt(2.0*_Ek/mHe);
    u0[9] = 0.0;
    u0[10] = 0.0;
    u0[11] = 0.0;

    double hmax = 1.0/(-u0[8]);
    double hmin = hmax/1.0e12;

    ODErkck().solve(*this, 0.0, 1.0e12, u0, 1.0e-8, hmin, hmax, _soln);

    final();
  }

  void d(double _t, const ODErkck::Data& _u, ODErkck::Data& _du) {
    /*** u = (x,y,z,X,Y,Z,dx,dy,dz,dX,dY,dZ) ***/
    double xyzHe[3] = {_u[0],_u[1],_u[2]};
    double xyzC[3] = {_u[3],_u[4],_u[5]};
    double xyz[3], r;
    double V, dV;
    int j, k, natoms = mol.natoms();

    mc.dEnergyMoveAtom(atomID,xyzC);
    for(k = 0; k < 6; k++) _du[k] = _u[k+6];
    for(k = 6; k < 12; k++) _du[k] = 0.0;
    for(k = 0; k < natoms; k++) {
      Util::vsub(xyzHe, mol.atom(k)->x, xyz, 3);
      r = Util::norm2(xyz,3);
      VLJHeC(r,V,dV);
      for(j = 0; j < 3; j++) _du[j+6] += -dV/r*xyz[j]/mHe;
      if(k == atomID)
	for(j = 0; j < 3; j++) _du[j+9] += dV/r*xyz[j]/mCarbon;
    }
    
    const double step = 1.0e-6;
    double new_x[3], xbak;
    double e1, e2, force;
    for(j = 0; j < 3; j++) new_x[j] = mol.atom(atomID)->x[j];
    for(j = 0; j < 3; j++) {
      xbak = new_x[j];
      new_x[j] = xbak-step;
      e1 = mc.dEnergyMoveAtom(atomID,new_x);
      mc.unMoveAtom();
      new_x[j] = xbak+step;
      e2 = mc.dEnergyMoveAtom(atomID,new_x);
      mc.unMoveAtom();
      force = -(e2-e1)/2.0/step;
      _du[j+9] += force/mCarbon;
      new_x[j] = xbak;
    }
  }

  bool stop(double _t, const ODErkck::Data& _u) {
    double r, xyz[3], xyzHe[3] = {_u[0],_u[1],_u[2]};
    bool stopFlag = true;
    int k, natoms = mol.natoms();
    for(k = 0; k < natoms; k++) {
      Util::vsub(xyzHe, mol.atom(k)->x, xyz, 3);
      r = Util::norm2(xyz,3);
      stopFlag = (stopFlag&&(r>2.0*dFar));
    }
    return stopFlag;
  }

private:
  AMod::Molecule mol;
  int atomID;
  double mHe, mCarbon, dFar;
  MC::MC_LCBOPII mc;

  void init(const AMod::Molecule& _mol, int _atomID) {
    mol = _mol;
    atomID = _atomID;
    mHe = 414.841;
    mCarbon = 1244.82;
    dFar = 10.0;
    mc.attach(mol);
    mc.init();
  }

  void final() { mc.final(); }

  double VLJ(double _e, double _u) {
    return 4.0*_e*(pow(_u,12)-pow(_u,6));
  }

  double dVLJ(double _e, double _u) {
    return 4.0*_e*(12*pow(_u,11)-6*pow(_u,5));
  }

  void VLJHeC(double _r, double& _V, double& _dV) {
    double e = 0.0014;
    double s = 2.74;
    double c = 2.5;
    
    if(_r > c*s) { 
      _V = 0.0; 
      _dV = 0.0;
    }
    else {
      _V = VLJ(e,s/_r)-VLJ(e,1.0/c);
      _dV = dVLJ(e,s/_r)*(-s/(_r*_r));
    }
  }
};

int main(int argc, char* argv[]) {
  if(argc != 6) {
    cout << "Usage: testScattHeC <file_mol> <atomID> <Ek> <s> <q>" << endl;
    return 1;
  }
  
  ifstream finMol;
  int atomID;
  double Ek, s, q;

  finMol.open(argv[1]);
  if(!finMol.good()) {
    cout << "Error: cannot open file " << argv[1] << endl;
    return 1;
  }
  if(!Util::readWord(argv[2],atomID)) {
    cout << "Error: cannot read <atomID> " << argv[2] << endl;
    return 2;
  }
  if(!Util::readWord(argv[3],Ek)) {
    cout << "Error: cannot read <Ek> " << argv[3] << endl;
    return 3;
  }
  if(!Util::readWord(argv[4],s)) {
    cout << "Error: cannot read <s> " << argv[4] << endl;
    return 4;
  }
  if(!Util::readWord(argv[5],q)) {
    cout << "Error: cannot read <q> " << argv[5] << endl;
    return 5;
  }

  AMod::Molecule mol;

  finMol >> mol;
  finMol.close();
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 6;
  }
  if(atomID < 0 || atomID >= mol.natoms()) {
    cout << "Error: illegal atomID " << atomID << endl;
    return 7;
  }

  HeCScatt::ODErkck::Soln soln;

  HeCScatt().orbit(mol,atomID,Ek,s,q,soln);

  AMod::MolO::setFormat(cout,12);
  for(int k = 0, kmax = soln.t.size(); k < kmax; k++) {
    cout << soln.t[k] << " ";
    for(int j = 0; j < 12; j++) cout << soln.u[k][j] << " ";
    cout << endl;
  }

  return 0;
}
