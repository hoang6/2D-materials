#include "../MC/MC_REBO.h"
#include "../MC/MC_LCBOPII.h"
#include "../Util/FuncDFunc.h"
#include "../Util/ConjGrad.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
using namespace std;

#define MC_TYPE MC_LCBOPII

class GPHMaker {
public:
  GPHMaker(int _ncell1D);
  ~GPHMaker();
  double make(double _a);
  double func(double* _pa);
  void dfunc(double* _pa, double* _pd);

  Util::Accessors<AMod::Molecule> molecule;

private:
  int ncell1D;
};

int main(int argc, char* argv[]) {
  if(argc != 5) {
    cout << "Usage: makeGPHGroundState <a0> <ncell1D> <auto_a/fix_a> <file_mol>" << endl;
    return 1;
  }

  /********** READ ARGUMENTS **********/
  ofstream mfout;
  istringstream iss;
  double a0;
  int ncell1D;
  int opt_a;
  
  AMod::MolO::setFormat(cout, 14);
  AMod::MolO::setFormat(mfout, 14);
  iss.str(string(argv[1])+" "+argv[2]+" "+argv[3]);
  iss >> a0 >> ncell1D >> opt_a;
  if(iss.fail()) {
    cout << "Error: illegal arguments" << endl;
    return 2;
  }
  mfout.open(argv[4]);
  if(!mfout.good()) {
    cout << "Error: cannot open " << argv[1] << endl;
    return 3;
  }

  /********** MAKE GPH **********/
  GPHMaker maker(ncell1D);
  if(opt_a == 0) {//auto_a
    int iterNum;
    double minE;
    double a = a0;
    Util::FuncDFunc<GPHMaker, double> fdf;
    fdf.init(maker, 1, &GPHMaker::func, &GPHMaker::dfunc, NULL, &cout);
    Util::conjGrad(&a, 1, 1.0e-9, 1000, fdf, iterNum, minE);
    maker.molecule().io().dumpTxt(mfout);
    cout << "Optimized a = " << a << endl;
  }
  else if(opt_a == 1) {//fix_a
    maker.make(a0);
    maker.molecule().io().dumpTxt(mfout);
  }
  cout << "natoms = " << maker.molecule().natoms() << endl;

  mfout.close();
  return 0;
}

/************************************************************/
GPHMaker::GPHMaker(int _ncell1D): ncell1D(_ncell1D) { molecule().reset(); }

GPHMaker::~GPHMaker() {}

double GPHMaker::make(double _a) {
  /********** MAKE GPH **********/
  int i, j, k, count;
  AMod::Atom::Data *pdataA, *pdataB;
  AMod::Molecule& _mol = molecule();
  double cell_a[] = {_a, 0.0, 0.0};
  double cell_b[] = {_a*0.5, _a*sqrt(3.0)/2.0, 0.0};
  double disp[] = {_a*0.5, _a/(2.0*sqrt(3.0)), 0.0};
  int _n = ncell1D;
  _mol.clear();
  _mol.axes().resize(2);
  _mol.atomsData().resize(_n*_n*2);
  for(k = 0; k < 3; k++) {
    _mol.axis(0).x[k] = _n*cell_a[k];
    _mol.axis(1).x[k] = _n*cell_b[k];
  }
  count = 0;
  for(i = 0; i < _n; i++) {
    for(j = 0; j < _n; j++) {
      pdataA = &_mol.atomsData()[count++];
      pdataB = &_mol.atomsData()[count++];
      pdataA->type = AMod::PTE::Carbon;
      pdataB->type = AMod::PTE::Carbon;
      for(k = 0; k < 3; k++) {
	pdataA->x[k] = i*cell_a[k]+j*cell_b[k];
	pdataB->x[k] = pdataA->x[k]+disp[k];
      }
    }//end for(int j...
  }//end for(int i...
  _mol.topo().sync();
  _mol.topo().update(AMod::CASE_AXES);
  _mol.topo().update(AMod::CASE_LONELY);

  /********** COMPUTE ENERGY **********/
  MC::MC_TYPE mc;
  mc.attach(_mol);
  mc.init();
  mc.final();
  return _mol.potential();
}

double GPHMaker::func(double* _pa) {
  return make(*_pa)/molecule().natoms();
}

void GPHMaker::dfunc(double* _pa, double* _pd) {
  const double step = 1.0e-6;
  double a, e1, e2;
  a = *_pa+step;
  e1 = this->func(&a);
  a = *_pa-step;
  e2 = this->func(&a);
  *_pd = (e1-e2)/2.0/step;
}


