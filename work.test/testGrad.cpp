#include "../AMod/AMod.h"
#include "../Util/util.h"
#include "../Util/RNG.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {
  if(argc != 2) {
    cout << "Usage: testGrad <file_mol>" << endl;
    return 1;
  }

  /********** Parse Arguments **********/  
  AMod::Molecule mol;
  ifstream fin;

  fin.open(argv[1]);
  if(!fin.good()) {
    cout << "Error: cannot open file: " << argv[1] << endl;
    return 2;
  }
  fin >> mol;
  if(mol.natoms() == 0) {
    cout << "Error: no molecule is read" << endl;
    return 3;
  }
  fin.close();
  mol.topo().mode(AMod::MODE_FREEZE);
  mol.topo().update(AMod::CASE_CHANGE_MODE);
  
  /********** compute gradient **********/
  AMod::MolOptim mopt(mol);
  const double step = 1.0e-6;

  mopt.init(AMod::REBO);
  int k, n = mopt.nDOF();
  double xbak, e1, e2;
  double* x = new double[n];
  double* g = new double[n];
  double* gref = new double[n];
  mopt.get_x(x);
  mopt.solver(x, gref);
  for(k = 0; k < n; k++) {
    if(k < 9 && mol.axis(k/3).fixed[k%3]) continue;
    if(k >= 9 && mol.atom((k-9)/3)->fixed[(k-9)%3]) continue;
    xbak = x[k];
    x[k] = xbak-step;
    e1 = mopt.solver(x, NULL);
    x[k] = xbak+step;
    e2 = mopt.solver(x, NULL);
    g[k] = (e2-e1)/2.0/step;
    x[k] = xbak;
  }
  mopt.final();

  /********** check gradient **********/
  double err, mse;
  mse = 0.0;
  AMod::MolO::setFormat(cout);
  for(k = 0; k < n; k++) {
    err = Util::abs(g[k]-gref[k]);
    cout << g[k] << " " << gref[k] << " " << (err>1.0e-3?"*":"") << endl;
    mse += err;
  }
  mse = sqrt(mse/n);

  cout << endl;
  cout << "mse = " << mse << endl;

  delete [] gref;
  delete [] g;
  delete [] x;

  return 0;
}
