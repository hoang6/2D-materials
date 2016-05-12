#include "../AMod/AMod.h"
#include "../MC/MC_REBO.h"
#include "../MC/MC_LCBOPII.h"
#include "../Util/util.h"
#include "../NR/jacobi.h"
#include "../NR/eigsrt.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <ctime>
using namespace std;

#define MC_TYPE MC_LCBOPII

typedef complex<double> DComplex;

void get3DHookMat(const AMod::Molecule& _mol, 
		  const double& _ddisp,
		  const double& _dstep,
		  vector<double>& _hookMat);

void getGPHDynMat(const AMod::Molecule& _mol, 
		  const vector<double>& _hookMat,
		  const double _qvec[2],
		  vector<DComplex>& _dynMat);

void eigGPHDynMat(const vector<DComplex>& _dynMat,
		  vector<double>& _d);

void dumpGPHDynMat(const AMod::Molecule& _mol, 
		   const vector<double>& _hookMat,
		   const double _ptStart[2],
		   const double _ptEnd[2],
		   const int& _npts,
		   bool _verbose,
		   ostream& _fout);

int main(int argc, char* argv[]) {
  if(argc != 8) {
    cout << "Usage: getGPHDynMat <file_mol> <delta_disp> <delta_step> <high_symm/mesh> <nq1D[0]> <nq1D[1]> <verbose>" << endl;
    cout << "Note: (delta_disp, delta_step) should be carefully chosen, e.g. (1.0e-5, 1.0e-6)" << endl;
    return 1;
  }
  
  /********** READ ARGUMENTS **********/
  AMod::Molecule mol;
  istringstream iss;
  double delta_disp;
  double delta_step;
  int opt_sample;
  int nq1D[2];
  bool verbose;
  ifstream mfin;

  AMod::MolO::setFormat(cout, 14);
  iss.str(string(argv[2])+" "+argv[3]+" "+argv[4]+" "+
	  argv[5]+" "+argv[6]+" "+argv[7]);
  iss >> delta_disp >> delta_step >> opt_sample 
      >> nq1D[0] >> nq1D[1] >> verbose;
  if(iss.fail()) {
    cout << "Error: illegal arguments" << endl;
    return 2;
  }
  mfin.open(argv[1]);
  if(!mfin.good()) {
    cout << "Error: cannot open file: " << argv[1] << endl;
    return 3;
  }
  mfin >> mol;
  if(mfin.fail() || mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 4;
  }
  mfin.close();

  /********** GET DYNAMIC MATRIX **********/
  double a;
  double qG[2], qM[2], qK[2];
  vector<double> hookMat;
  a = sqrt(3.0*Util::vvdist2(mol.atom(0)->x, mol.atom(1)->x, 3));
  qG[0] = 0.0;                  qG[1] = 0.0;
  qM[0] = Util::PI/a;           qM[1] = -Util::PI/(sqrt(3.0)*a);
  qK[0] = Util::PI*4.0/(3.0*a); qK[1] = 0.0;
  get3DHookMat(mol, delta_disp, delta_step, hookMat);
  /* dump hookMat */
  if(verbose) {
    for(int k = 0; k < int(hookMat.size()); k++) 
      cout << hookMat[k] << " ";
    cout << endl << endl;
  }
  if(opt_sample == 0) {//high symm
    int nq1DGM, nq1DMK, nq1DKG;
    nq1DGM = int(nq1D[0]/2.0);
    nq1DMK = int(nq1D[0]/(2.0*sqrt(3.0)));
    nq1DKG = int(nq1D[0]/(sqrt(3.0)));
    //Gamma-->M    
    dumpGPHDynMat(mol, hookMat, qG, qM, nq1DGM, verbose, cout);
    cout << endl;
    //M-->K
    dumpGPHDynMat(mol, hookMat, qM, qK, nq1DMK, verbose, cout);
    cout << endl;
    //K-->Gamma
    dumpGPHDynMat(mol, hookMat, qK, qG, nq1DKG, verbose, cout);
    cout << endl;
  }
  else if(opt_sample == 1) {//mesh
    double qRadius = 4.0*Util::PI/(sqrt(3.0)*a);
    double qX[2] = {sqrt(3.0)/2.0*qRadius, -0.5*qRadius};
    double qY[2] = {0.0, qRadius};
    double qA[2], qB[2];
    for(int k = 0; k < nq1D[0]; k++) {
      for(int j = 0; j < 2; j++) {
	qA[j] = k*qX[j]/nq1D[0];
	qB[j] = qA[j]+qY[j];
      }
      //A-->B
      dumpGPHDynMat(mol, hookMat, qA, qB, nq1D[1], verbose, cout);
    }//end for(int k...
    cout << endl;
  }

  return 0;
}

void get3DHookMat(const AMod::Molecule& _mol, 
		  const double& _ddisp,
		  const double& _dstep,
		  vector<double>& _hookMat) {
  AMod::Molecule mol;
  MC::MC_TYPE mc;
  double dE, f, new_x[3];
  int i, j, k, n, ai, aj;
  
  mol = _mol;
  mc.attach(mol);
  mc.init();  
  n = 3*mol.natoms();
  _hookMat.assign(6*n, 0.0);
  for(k = -1; k <= 1; k+=2) {//replica
    for(i = 0; i < 6; i++) {
      mc.save();
      ai = i/3;
      Util::vcopy(mol.atom(ai)->x, new_x, 3);
      new_x[i%3] += k*_ddisp;
      mc.dEnergyMoveAtom(ai, new_x);
      for(j = 0; j < n; j++) {
	aj = j/3;
	Util::vcopy(mol.atom(aj)->x, new_x, 3);
	new_x[j%3] += _dstep;
	dE = mc.dEnergyMoveAtom(aj, new_x);
	mc.unMoveAtom();
	new_x[j%3] -= 2.0*_dstep;
	dE -= mc.dEnergyMoveAtom(aj, new_x);
	mc.unMoveAtom();
	f = -dE/(2.0*_dstep);
	_hookMat[i+j*6] += -f/(k*_ddisp)/2.0;
      }//end for(j...
      mc.restore();
    }//end for(i...
  }//end for(k...
  mc.final();
}

void getGPHDynMat(const AMod::Molecule& _mol, 
		  const vector<double>& _hookMat, 
		  const double _qvec[2],
		  vector<DComplex>& _dynMat) {
  int i, j, jp, n;
  double dR[3];
  double phase;

  n = 3*_mol.natoms();
  _dynMat.assign(6*6, DComplex(0.0,0.0));
  for(i = 0; i < 6; i++) {
    for(j = 0; j < 6; j++) {
      DComplex& entry = _dynMat[i+j*6];
      for(jp = j; jp < n; jp += 6) {
	_mol.topo().displacement(0, int(jp/6)*2, dR);
	phase = -(_qvec[0]*dR[0]+_qvec[1]*dR[1]);
	entry += _hookMat[i+jp*6]*DComplex(cos(phase),sin(phase));
      }//end for(jp...
    }//end for(j...
  }//end for(i...
}

void eigGPHDynMat(const vector<DComplex>& _dynMat,
		  vector<double>& _d) {
  int n = 6;
  NR::Matrix<double> a(1,2*n,1,2*n);
  NR::Matrix<double> v(1,2*n,1,2*n);
  NR::Vector<double> d(1,2*n);
  double** pa = a.pointer();
  double** pv = v.pointer();
  double* pd = d.pointer();
  int i, j, nrot;
  for(i = 1; i <= 2*n; i++) {
    for(j = i; j <= 2*n; j++) {
      if(i <= n && j <= n) {
	if(i == j) 
	  pa[i][j] = _dynMat[(i-1)+(i-1)*n].real();
	else
	  pa[i][j] = (_dynMat[(i-1)+(j-1)*n].real()+
		      _dynMat[(j-1)+(i-1)*n].real())/2.0;
      }
      else if(i <= n && j > n) {
	if(i == j-n)
	  pa[i][j] = 0.0;
	else
	  pa[i][j] = (_dynMat[(j-n-1)+(i-1)*n].imag()-
		      _dynMat[(i-1)+(j-n-1)*n].imag())/2.0;
      }
      else {
	pa[i][j] = pa[i-n][j-n];
      }
    }//end for(j...
  }//end for(i...
  NR::jacobi(pa,2*n,pd,pv,&nrot);
  NR::eigsrt(pd,pv,2*n);
  _d.resize(n);
  for(i = 0; i < n; i++) _d[i] = pd[2*i+1];
}

void dumpGPHDynMat(const AMod::Molecule& _mol, 
		   const vector<double>& _hookMat,
		   const double _ptStart[2],
		   const double _ptEnd[2],
		   const int& _nq1D,
		   bool _verbose,
		   ostream& _fout) {
  int j, k, n = 6;
  double qvec[2];
  vector<DComplex> dynMat;
  vector<double> d;
  if(_nq1D <= 0) {
    _fout << _ptStart[0] << " " << _ptStart[1] << " ";
    if(_verbose)
      for(j = 0; j < 2*n*n+n; j++) _fout << 0.0 << " ";
    else
      for(j = 0; j < n; j++) _fout << 0.0 << " ";
    _fout << endl;
    return;
  }
  for(k = 0; k < _nq1D; k++) {
    for(j = 0; j < 2; j++) 
      qvec[j] = _ptStart[j]+k*(_ptEnd[j]-_ptStart[j])/double(_nq1D);
    getGPHDynMat(_mol, _hookMat, qvec, dynMat);
    eigGPHDynMat(dynMat, d);
    _fout << qvec[0] << " " << qvec[1] << " ";
    if(_verbose) {
      for(j = 0; j < n*n; j++) 
	_fout << dynMat[j].real() << " " << dynMat[j].imag() << " ";
    }
    for(j = 0; j < n; j++) _fout << d[j] << " ";
    _fout << endl;
  }//end for(k...
}
